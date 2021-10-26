using Base.Order



"""
    serializable(order::Ordering, x)

Wheather [`serialize`](@ref) and [`deserialize`](@ref) are implemented for the
provided arguments. In the future, there may be an option for serializable but not
deserializable
"""
serializable(order::Ordering, T::Type) = serialized_type(order, T) != Nothing

"""
    serialized_type(order::Ordering, T::Type)

The return type of [`serialize(o, x::T)`](@ref serialize) if `serialize` and
[`deserialize`](@ref) are implemented for the provided arguments. If not, implemented,
return Nothing.
"""
serialized_type(order::Ordering, x::Type) = typeof(serialize(order, zero(x)))

"""
    serialize(order::Ordering, x)::Unsigned

Map `x` to an un unsigned integer, maintaining sort order.

The map should be reversible with [`deserialize`](@ref), so `!lt(order, a, b)` must be
total for `a, b <: typeof(x)`. Satisfies
`lt(order, a, b) === (serialize(order, a) < serialize(order, b))`
and `x === deserialize(T, order, serialize(order, x::T))`

EXCEPTION: scalar serialize fails on NaN edge cases. This is handled in vector serialize.

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
serialized_type(::ForwardOrdering, T::Type{<:Union{Unsigned, Signed}}) = isbitstype(T) ? unsigned(T) : Nothing

# Floats
for (float, int) in ((Float32, Int32), (Float64, Int64))
    @eval function serialize(::ForwardOrdering, x::$float)
        y = reinterpret($int, x)
        unsigned(y < 0 ? ~y : xor(y, typemin(y))) - ~reinterpret(unsigned($int), $float(-Inf))
    end
    @eval function deserialize(T::Type{$float}, ::ForwardOrdering, u::unsigned($int))
        y = reinterpret($int, u + ~reinterpret(unsigned($int), $float(-Inf)))
        reinterpret(T, y < 0 ? xor(y, typemin(y)) : ~y)
    end
end

# Booleans
serialize(::ForwardOrdering, x::Bool) = UInt8(x)
deserialize(::Type{Bool}, ::ForwardOrdering, u::UInt8) = Bool(u)
serialized_type(::ForwardOrdering, ::Type{Bool}) = UInt8

# Chars
serialize(::ForwardOrdering, x::Char) = UInt32(x)
deserialize(::Type{Char}, ::ForwardOrdering, u::UInt32) = Char(u)
serialized_type(::ForwardOrdering, ::Type{Char}) = UInt32


### Reverse orderings
serialize(rev::ReverseOrdering, x) = ~serialize(rev.fwd, x)
deserialize(T::Type, rev::ReverseOrdering, u::Unsigned) = deserialize(T, rev.fwd, ~u)
serialized_type(order::ReverseOrdering, T::Type) where F = serialized_type(order.fwd, T)


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
serialize!(xs::AbstractVector, lo::Integer, hi::Integer, order::Ordering) =
    serialize!(reinterpret(serialized_type(order, eltype(xs))), xs, lo, hi, order)
serialize(xs::AbstractVector, lo::Integer, hi::Integer, order::Ordering) =
    serialize!(OffsetVector{serialized_type(order, eltype(xs))}(undef, axes(xs)), xs, lo, hi, order)

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
