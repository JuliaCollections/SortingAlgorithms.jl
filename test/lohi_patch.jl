radix_sort!(ps::AbstractVector{U}, ts::AbstractVector{U}, ::Unsigned, bits::Unsigned, chunk_size::Unsigned) where {U <: Unsigned} =
    radix_sort!(ps, ts, firstindex(ps), lastindex(ps), bits, chunk_size)

serialize!(order::Ordering, xs::AbstractVector) =
    serialize!(xs, firstindex(xs), lastindex(xs), order)
serialize(order::Ordering, xs::AbstractVector) =
    serialize(xs, firstindex(xs), lastindex(xs), order)

deserialize!(xs::AbstractVector, order::Ordering, us::AbstractVector{<:Unsigned}, compression=nothing) =
    deserialize!(xs, us, firstindex(xs), lastindex(xs), order, compression)
