using Base.Order



"""
    Serializable(order::Ordering, T::Type)

Return `Some(typeof(serialize(o, ::T)))` if [`serialize`](@ref) and
[`deserialize`](@ref) are implemented.

If either is not implemented, return nothing.

EXCEPTION: For containters, return the serialized element type.
"""
Serializable(order::Ordering, T::Type) = nothing

"""
    serialize(order::Ordering, x)::Unsigned

Map `x` to an un unsigned integer, maintaining sort order.

The map should be reversible with [`deserialize`](@ref), so `!lt(order, a, b)` must be
total for `a, b <: typeof(x)`. Satisfies
`lt(order, a, b) === (serialize(order, a) < serialize(order, b))`
and `x === deserialize(T, order, serialize(order, x::T))`

NOTE: scalar serialize fails on NaN edge cases. This is handled in vector serialize.
Consequently `Serialize(Float64) === nothing`.

See also: [`serializable`](@ref) [`deserialize`](@ref)
"""
function serialize end

"""
    deserialize(T::Type, order::Ordering, u::Unsigned)

Reconstruct the unique value `x::T` that serializes to `u`. Satisfies
`x === deserialize(T, order, serialize(order, x::T))` for all `x <: T`.

See also: [`serialize`](@ref) [`serializable`](@ref)
"""
function deserialize end


### Primitive Types

# Integers
serialize(::ForwardOrdering, x::Unsigned) = x
deserialize(::Type{T}, ::ForwardOrdering, u::T) where T <: Unsigned = u
serialize(::ForwardOrdering, x::Signed) = unsigned(xor(x, typemin(x)))
deserialize(::Type{T}, ::ForwardOrdering, u::Unsigned) where T <: Signed = xor(signed(u), typemin(T))
Serializable(::ForwardOrdering, T::Type{<:Union{Unsigned, Signed}}) = isbitstype(T) ? Some(unsigned(T)) : nothing

# Floats are not Serializable under regular orderings because they fail on NaN edge cases.
for (float, int) in ((Float16, Int16), (Float32, Int32), (Float64, Int64))
    @eval function serialize(::Base.Sort.Float.Right, x::$float)
        y = reinterpret($int, x)
        unsigned(y < 0 ? ~y : xor(y, typemin(y))) - ~reinterpret(unsigned($int), $float(-Inf))
    end
    @eval function deserialize(T::Type{$float}, ::Base.Sort.Float.Right, u::unsigned($int))
        y = reinterpret($int, u + ~reinterpret(unsigned($int), $float(-Inf)))
        reinterpret(T, y < 0 ? xor(y, typemin(y)) : ~y)
    end
    @eval Serializable(::Base.Sort.Float.Right, ::Type{$float}) = unsigned($int)

    #This repetitive code is necessary because we have Base.Sort.Float.Left rather than
    #Base.Order.ReverseOrdering{Base.Sort.Float.Right}.
    @eval function serialize(::Base.Sort.Float.Left, x::$float)
        y = reinterpret($int, x)
        unsigned(y < 0 ? ~y : xor(y, typemin(y))) - ~reinterpret(unsigned($int), $float(-Inf))
    end
    @eval function deserialize(T::Type{$float}, ::Base.Sort.Float.Left, u::unsigned($int))
        y = reinterpret($int, u + ~reinterpret(unsigned($int), $float(-Inf)))
        reinterpret(T, y < 0 ? xor(y, typemin(y)) : ~y)
    end
    @eval Serializable(::Base.Sort.Float.Left, ::Type{$float}) = unsigned($int)
end

# Booleans
serialize(::ForwardOrdering, x::Bool) = UInt8(x)
deserialize(::Type{Bool}, ::ForwardOrdering, u::UInt8) = Bool(u)
Serializable(::ForwardOrdering, ::Type{Bool}) = UInt8

# Chars
serialize(::ForwardOrdering, x::Char) = reinterpret(UInt32, x)
deserialize(::Type{Char}, ::ForwardOrdering, u::UInt32) = reinterpret(Char, u)
Serializable(::ForwardOrdering, ::Type{Char}) = UInt32


### Reverse orderings
serialize(rev::ReverseOrdering, x) = ~serialize(rev.fwd, x)
deserialize(T::Type, rev::ReverseOrdering, u::Unsigned) = deserialize(T, rev.fwd, ~u)
Serializable(order::ReverseOrdering, T::Type) = Serializable(order.fwd, T)


### Vectors
"""docstring"""
function serialize!(us::AbstractVector{<:Unsigned}, xs::AbstractVector, lo::Integer, hi::Integer, order::Ordering)
    us[lo] = mn = mx = serialize(order, xs[lo])
    i = lo # rename lo -> i for clarity only
    while i < hi
        i += 1

        u = us[i] = serialize(order, xs[i])

        if u > mx
            mx = u
        elseif u < mn
            mn = u
        end
    end
    us, mn, mx
end
serialize!(T::Type{<:Unsigned}, xs::AbstractVector, lo::Integer, hi::Integer, order::Ordering) =
    serialize!(reinterpret(T, xs), xs, lo, hi, order)
serialize(T::Type{<:Unsigned}, xs::AbstractVector, lo::Integer, hi::Integer, order::Ordering) =
    serialize!(OffsetVector{T}(undef, axes(xs)), xs, lo, hi, order)

compress!(us::AbstractVector, ::Nothing) = us
compress!(us::AbstractVector{U}, min::U) where U <: Unsigned = us .-= min

"""docstring"""
function deserialize!(xs::AbstractVector, us::AbstractVector{<:Unsigned},
    lo::Integer, hi::Integer, order::Ordering, compression::Nothing)
    @inbounds for i in lo:hi
        xs[i] = deserialize(eltype(xs), order, us[i])
    end
    xs
end
function deserialize!(xs::AbstractVector, us::AbstractVector{U},
    lo::Integer, hi::Integer, order::Ordering, compression::U) where U <: Unsigned
    @inbounds for i in lo:hi
        xs[i] = deserialize(eltype(xs), order, us[i]+compression)
    end
    xs
end

Serializable(order::Base.Order.ReverseOrdering, T::Type{<:AbstractVector}) = Serializable(order.fwd, eltype(T))
Serializable(order::Ordering, T::Type{<:AbstractVector}) = Serializable(order, eltype(T))


### Notes

#TODO specialize vectorized serialization on floats follwing from fpsort!


#TODO Extend to other orderings and/or types



#= Notes on floating points
we could change ForwardOrdering => Union{Base.Sort.Float.Left, ForwardOrdering}
which would be consistant,
but we should never both serialize and use Base.Sort.Float.fpsort!
and indeed this whole thing should superseede Base.Sort.Float.fpsort!

TODO in Base/sort.jl if we want to comap with it:
- struct Right <: Ordering end
- right(::DirectOrdering) = Right()
+ right(::DirectOrdering) = Reverse(Left())
- lt(::Right, x::T, y::T) where {T<:Floats} = slt_int(x, y)
=#

#TODO types are too wide all over this file, especially of the form AbstractVector{<:Unsigned}.
#TODO ensure that reinterpret is fast in use


# Checking serializability, Three options:
# 1) Manual *current*
#   a) general definition of seializable as false
#   b) serialize is defined `function serialize end` without defaults
#   c) wherever serialize is specialized, so is serializable
#   +) performant, straightforward
#   -) more code, users need to specialize serializable & serialize
# 2) Run it
#   a) serialize is defined to regturn missing by default
#   b) serializable runs serialize concretely
#   +) concise, easy to extend serialize
#   -) *no support for eltype*, performance overhead
# 3) Muck with Base.return_types etc.
#   a) idk
#   +) easy for user to specialize serialize & supports eltype.
#   -) mucking around

#TODO integrate nan and missing filtration into vector serialization
#(requires passing a new hi back from serialize!)

# consider removing lo and hi in favor of view everywhere
# (note however: https://github.com/JuliaLang/julia/issues/39864)
# so for now we will pass lo and high everywhere :(
# this is exactly what views are for.

#TODO consider Float16, and deal with teh fact that scalar float serialization does not meet
# docstring specs.
