import Base: Forward, ForwardOrdering, Reverse, ReverseOrdering, Lexicographic, LexicographicOrdering, sortperm, Ordering, 
            setindex!, getindex, similar

# create bits types for easy loading of bytes of lengths up to 15
primitive type Bits24 24 end
primitive type Bits40 40 end
primitive type Bits48 48 end
primitive type Bits56 56 end
primitive type Bits72 72 end
primitive type Bits80 80 end
primitive type Bits88 88 end
primitive type Bits96 96 end
primitive type Bits104 104 end
primitive type Bits112 112 end
primitive type Bits120 120 end


function load_bits_with_padding(::Type{UInt128}, s::String, skipbytes = 0)::UInt128
    n = sizeof(s)

    remaining_bytes_to_load = min(sizeof(s), n - skipbytes)
    nbits_to_shift_away = 8(sizeof(UInt128) - remaining_bytes_to_load)

    ptrs  = pointer(s) + skipbytes

    # the below checks if the string is less than 8 bytes away from the page
    # boundary assuming the page size is 4kb
    # see https://discourse.julialang.org/t/is-there-a-way-to-check-how-far-away-a-pointer-is-from-a-page-boundary/8147/11?u=xiaodai
    if (UInt(ptrs) & 0xfff) > 0xff8
        if  remaining_bytes_to_load == 15
            return ntoh(Base.zext_int(UInt128, unsafe_load(Ptr{Bits120}(ptrs))))
        elseif  remaining_bytes_to_load == 14
            return ntoh(Base.zext_int(UInt128, unsafe_load(Ptr{Bits112}(ptrs))))
        elseif  remaining_bytes_to_load == 13
            return ntoh(Base.zext_int(UInt128, unsafe_load(Ptr{Bits104}(ptrs))))
        elseif  remaining_bytes_to_load == 12
            return ntoh(Base.zext_int(UInt128, unsafe_load(Ptr{Bits96}(ptrs))))
        elseif  remaining_bytes_to_load == 11
            return ntoh(Base.zext_int(UInt128, unsafe_load(Ptr{Bits88}(ptrs))))
        elseif  remaining_bytes_to_load == 10
            return ntoh(Base.zext_int(UInt128, unsafe_load(Ptr{Bits80}(ptrs))))
        elseif  remaining_bytes_to_load == 9
            return ntoh(Base.zext_int(UInt128, unsafe_load(Ptr{Bits72}(ptrs))))
        elseif  remaining_bytes_to_load == 8
            return ntoh(Base.zext_int(UInt128, unsafe_load(Ptr{UInt64}(ptrs))))
        elseif  remaining_bytes_to_load == 7
            return ntoh(Base.zext_int(UInt128, unsafe_load(Ptr{Bits56}(ptrs))))
        elseif  remaining_bytes_to_load == 6
            return ntoh(Base.zext_int(UInt128, unsafe_load(Ptr{Bits48}(ptrs))))
        elseif  remaining_bytes_to_load == 5
            return ntoh(Base.zext_int(UInt128, unsafe_load(Ptr{Bits40}(ptrs))))
        elseif  remaining_bytes_to_load == 4
            return ntoh(Base.zext_int(UInt128, unsafe_load(Ptr{UInt32}(ptrs))))
        elseif  remaining_bytes_to_load == 3
            return ntoh(Base.zext_int(UInt128, unsafe_load(Ptr{Bits24}(ptrs))))
        elseif  remaining_bytes_to_load == 2
            return ntoh(Base.zext_int(UInt128, unsafe_load(Ptr{UInt16}(ptrs))))
        else
            return ntoh(Base.zext_int(UInt128, unsafe_load(Ptr{UInt8}(ptrs))))
        end
    else
        return ntoh(unsafe_load(Ptr{UInt128}(ptrs))) >> nbits_to_shift_away << nbits_to_shift_away
    end 
end

# these functions assumes that remaining_bytes_to_load > 0
function load_bits_with_padding(::Type{UInt64}, s::String, skipbytes = 0)::UInt64
    n = sizeof(s)

    remaining_bytes_to_load = min(sizeof(s), n - skipbytes)
    nbits_to_shift_away = 8(sizeof(UInt64) - remaining_bytes_to_load)

    ptrs  = pointer(s) + skipbytes
    # the below checks if the string is less than 8 bytes away from the page
    # boundary assuming the page size is 4kb
    # see https://discourse.julialang.org/t/is-there-a-way-to-check-how-far-away-a-pointer-is-from-a-page-boundary/8147/11?u=xiaodai
    if (UInt(ptrs) & 0xfff) > 0xff8
        if  remaining_bytes_to_load == 7
            return ntoh(Base.zext_int(UInt, unsafe_load(Ptr{Bits56}(ptrs))))
        elseif  remaining_bytes_to_load == 6
            return ntoh(Base.zext_int(UInt, unsafe_load(Ptr{Bits48}(ptrs))))
        elseif  remaining_bytes_to_load == 5
            return ntoh(Base.zext_int(UInt, unsafe_load(Ptr{Bits40}(ptrs))))
        elseif  remaining_bytes_to_load == 4
            return ntoh(Base.zext_int(UInt, unsafe_load(Ptr{UInt}(ptrs))))
        elseif  remaining_bytes_to_load == 3
            return ntoh(Base.zext_int(UInt, unsafe_load(Ptr{Bits24}(ptrs))))
        elseif  remaining_bytes_to_load == 2
            return ntoh(Base.zext_int(UInt, unsafe_load(Ptr{UInt16}(ptrs))))
        else
            return ntoh(Base.zext_int(UInt, unsafe_load(Ptr{UInt8}(ptrs))))
        end
    else
        return ntoh(unsafe_load(Ptr{UInt64}(ptrs))) >> nbits_to_shift_away << nbits_to_shift_away
    end 
end

function load_bits_with_padding(::Type{UInt32}, s::String, skipbytes = 0)::UInt32
    n = sizeof(s)

    remaining_bytes_to_load = min(sizeof(s), n - skipbytes)
    nbits_to_shift_away = 8(sizeof(UInt32) - remaining_bytes_to_load)

    ptrs  = pointer(s) + skipbytes
    # the below checks if the string is less than 8 bytes away from the page
    # boundary assuming the page size is 4kb
    # see https://discourse.julialang.org/t/is-there-a-way-to-check-how-far-away-a-pointer-is-from-a-page-boundary/8147/11?u=xiaodai
    if (UInt(ptrs) & 0xfff) > 0xff8
        if  remaining_bytes_to_load == 3
            return ntoh(Base.zext_int(UInt, unsafe_load(Ptr{Bits24}(ptrs))))
        elseif  remaining_bytes_to_load == 2
            return ntoh(Base.zext_int(UInt, unsafe_load(Ptr{UInt16}(ptrs))))
        else
            return ntoh(Base.zext_int(UInt, unsafe_load(Ptr{UInt8}(ptrs))))
        end
    else
        return ntoh(unsafe_load(Ptr{UInt32}(ptrs))) >> nbits_to_shift_away << nbits_to_shift_away
    end 
end

"""
    load_bits([type,] s, skipbytes)

Load the underlying bits of a string `s` into a `type` of the user's choosing.
The default is `UInt`, so on a 64 bit machine it loads 64 bits (8 bytes) at a time.
If the `String` is shorter than 8 bytes then it's padded with 0.

- `type`:       any bits type that has `>>`, `<<`, and `&` operations defined
- `s`:          a `String`
- `skipbytes`:  how many bytes to skip e.g. load_bits("abc", 1) will load "bc" as bits
"""
# Some part of the return result should be padded with 0s.
# To prevent any possibility of segfault we load the bits using
# successively smaller types
# it is assumed that the type you are trying to load into needs padding
# i.e. `remaining_bytes_to_load > 0`
function load_bits(::Type{T}, s::String, skipbytes = 0)::T where T <: Unsigned
    n = sizeof(s)
    if n < skipbytes
        res = zero(T)
    elseif n - skipbytes >= sizeof(T)
        res = ntoh(unsafe_load(Ptr{T}(pointer(s, skipbytes+1))))
    else
        res = load_bits_with_padding(T, s, skipbytes)
    end
    return res
end

# Radix sort for strings
function sort!(svec::AbstractVector{String}, ::StringRadixSortAlg, o::Perm)
    sort!(svec, 1, length(svec), StringRadixSort, o)
end

function sort!(svec::AbstractVector, lo::Int, hi::Int, ::StringRadixSortAlg, o::O) where O <: Union{ForwardOrdering, ReverseOrdering, Perm}
    if isa(o, Perm)
        if eltype(o.data) != String
            throw(ArgumentError("Cannot use StringRadixSort on type $(eltype(o.data))"))
        end
        return sortperm_radixsort(o.data, order = o.order)
    else
        sort!(svec, lo, hi, StringRadixSort, o)
    end
end

function sort!(svec::AbstractVector{String}, lo::Int, hi::Int, ::StringRadixSortAlg, o::O) where O <: Union{ForwardOrdering, ReverseOrdering}
    #  Input checking
    # this should never be the case where `o isa Perm` but the `svec`` is a string
    # if isa(o, Perm)
    #     throw(ArgumentError("Cannot use StringRadixSort on type $(eltype(o.data))"))
        # if eltype(o.data) != String
        #     throw(ArgumentError("Cannot use StringRadixSort on type $(eltype(o.data))"))
        # end
        # o = o.order
        # svec = o.data
        # return
    # end
    
    if lo >= hi;  return svec;  end
    # the length subarray to sort
    l = hi - lo + 1

    # find the maximum string length    
    lens = maximum(sizeof, svec)
    skipbytes = lens
    if lens > 0
        while lens > 8
            skipbytes = max(0, skipbytes - 16)
            bits128 = zeros(UInt128, l)
            if o == Reverse
                bits128[lo:hi] .= .~load_bits.(UInt128, @view(svec[lo:hi]), skipbytes)
                sorttwo!(bits128, svec, lo, hi)
            else
                bits128[lo:hi] .= load_bits.(UInt128, @view(svec[lo:hi]), skipbytes)
                sorttwo!(bits128, svec, lo, hi)
            end
            lens -= 16
        end
        if lens > 4
            skipbytes = max(0, skipbytes - 8)
            bits64 = zeros(UInt64, l)
            if o == Reverse
                bits64[lo:hi] .= .~load_bits.(UInt64, @view(svec[lo:hi]), skipbytes)
                sorttwo!(bits64, svec, lo, hi)
            else
                bits64[lo:hi] .= load_bits.(UInt64, @view(svec[lo:hi]), skipbytes)
                sorttwo!(bits64, svec, lo, hi)
            end
            lens -= 8
        end
        if lens > 0
            skipbytes = max(0, skipbytes - 4)
            bits32 = zeros(UInt32, l)
            if o == Reverse
                bits32[lo:hi] .= .~load_bits.(UInt32, @view(svec[lo:hi]), skipbytes)
                sorttwo!(bits32, svec, lo, hi)
            else
                bits32[lo:hi] .= .~load_bits.(UInt32, @view(svec[lo:hi]), skipbytes)
                sorttwo!(bits32, svec, lo, hi)
            end
            lens -= 4
        end
    end
    svec
end

"""
    sorttwo!(vs, index)

Sort both the `vs` and reorder `index` at the same. This allows for faster sortperm
for radix sort.
"""
function sorttwo!(vs::AbstractVector{T}, index, lo::Int = 1, hi::Int=length(vs)) where T
    # Input checking
    if lo >= hi;  return (vs, index);  end

    # Make sure we're sorting a bits type
    # T = Base.Order.ordtype(o, vs)
    if !isbits(T)
        error("Radix sort only sorts bits types (got $T)")
    end
    o = Forward

    # Init
    iters = ceil(Integer, sizeof(T)*8/RADIX_SIZE)
    # number of buckets in the counting step
    nbuckets = 2^RADIX_SIZE
    bin = zeros(UInt32, nbuckets, iters)
    if lo > 1;  bin[1,:] = lo-1;  end

    # Histogram for each element, radix
    for i = lo:hi
        v = uint_mapping(o, vs[i])
        for j = 1:iters
            idx = Int((v >> (j-1)*RADIX_SIZE) & RADIX_MASK) + 1
            @inbounds bin[idx,j] += 1
        end
    end

    # Sort!
    swaps = 0
    len = hi-lo+1

    index1 = similar(index)
    ts=similar(vs)
    for j = 1:iters
        # Unroll first data iteration, check for degenerate case
        v = uint_mapping(o, vs[hi])
        idx = Int((v >> (j-1)*RADIX_SIZE) & RADIX_MASK) + 1

        # are all values the same at this radix?
        if bin[idx,j] == len;  continue;  end

        # cbin = cumsum(bin[:,j])
        # tries to achieve the above one-liner with more efficiency
        cbin = zeros(UInt32, nbuckets)
        cbin[1] = bin[1,j]
        for i in 2:nbuckets
            cbin[i] = cbin[i-1] + bin[i,j]
        end

        ci = cbin[idx]
        ts[ci] = vs[hi]
        index1[ci] = index[hi]
        cbin[idx] -= 1

        # Finish the loop...
        @inbounds for i in hi-1:-1:lo
            v = uint_mapping(o, vs[i])
            idx = Int((v >> (j-1)*RADIX_SIZE) & RADIX_MASK) + 1
            ci = cbin[idx]
            ts[ci] = vs[i]
            index1[ci] = index[i]
            cbin[idx] -= 1
        end
        vs,ts = ts,vs
        index, index1 = index1, index
        swaps += 1
    end

    if isodd(swaps)
        vs,ts = ts,vs
        index, index1 = index1, index
        for i = lo:hi
            @inbounds vs[i] = ts[i]
            @inbounds index[i] = index1[i]
        end
    end
    (vs, index)
end

"""
    sortperm_radixsort(svec, rev = nothing, order = Forward)

To return a `String` vector using LSD radixsort
"""
function sortperm_radixsort(svec::AbstractVector{String}; rev::Union{Bool,Void}=nothing, order::Ordering=Forward)
    sortperm_radixsort!(copy(svec), rev = rev, order =order)
end

function sortperm_radixsort!(svec::AbstractVector{String}; rev::Union{Bool,Void}=nothing, order::Ordering=Forward)
    siv = StringIndexVector(svec, collect(1:length(svec)))

    # find the maximum string length
    lens = reduce((x,y) -> max(x,sizeof(y)),0, svec)
    skipbytes = lens
    while lens > 0
       if lens > 8
            skipbytes = max(0, skipbytes - 16)
            if order == Reverse
                sorttwo!(.~load_bits.(UInt128, siv.svec, skipbytes), siv)
            else
                sorttwo!(load_bits.(UInt128, siv.svec, skipbytes), siv)
            end
            lens -= 16
        elseif lens > 4
            skipbytes = max(0, skipbytes - 8)
            if order == Reverse
                sorttwo!(.~load_bits.(UInt64, siv.svec, skipbytes), siv)
            else
                sorttwo!(load_bits.(UInt64, siv.svec, skipbytes), siv)
            end
            lens -= 8
        else
            skipbytes = max(0, skipbytes - 4)
            if order == Reverse
                sorttwo!(.~load_bits.(UInt32, siv.svec, skipbytes), siv)
            else
                sorttwo!(load_bits.(UInt32, siv.svec, skipbytes), siv)
            end
            lens -= 4
        end
    end
    siv.index
end

"Simple data structure for carrying a string vector and its index; this allows
`sorttwo!` to sort the radix of the string vector and reorder the string and its
index at the same time opening the way for faster sort_perm"
struct StringIndexVector
    svec::Vector{String}
    index::Vector{Int}
end

function setindex!(siv::StringIndexVector, X::StringIndexVector, inds)
    siv.svec[inds] = X.svec
    siv.index[inds] = X.index
end

function setindex!(siv::StringIndexVector, X, inds)
    siv.svec[inds] = X[1]
    siv.index[inds] = X[2]
end

getindex(siv::StringIndexVector, inds::Integer) = siv.svec[inds], siv.index[inds]
getindex(siv::StringIndexVector, inds...) = StringIndexVector(siv.svec[inds...], siv.index[inds...])
similar(siv::StringIndexVector) = StringIndexVector(similar(siv.svec), similar(siv.index))