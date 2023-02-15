__precompile__()

module SortingAlgorithms

using DataStructures
using Base.Sort
using Base.Order
using Base: Cartesian

import Base.Sort: sort!
import DataStructures: heapify!, percolate_down!
import StaticArrays: MVector

export HeapSort, TimSort, RadixSort, CombSort, BranchyPatternDefeatingQuicksort, BranchlessPatternDefeatingQuicksort, BranchyPdqSort, BranchlessPdqSort

struct HeapSortAlg  <: Algorithm end
struct TimSortAlg   <: Algorithm end
struct RadixSortAlg <: Algorithm end
struct CombSortAlg  <: Algorithm end
abstract type PatternDefeatingQuicksortAlg  <: Algorithm end
struct BranchyPatternDefeatingQuicksortAlg  <: PatternDefeatingQuicksortAlg end
struct BranchlessPatternDefeatingQuicksortAlg  <: PatternDefeatingQuicksortAlg end

function maybe_optimize(x::Algorithm) 
    isdefined(Base.Sort, :InitialOptimizations) ? Base.Sort.InitialOptimizations(x) : x
end     
const HeapSort  = maybe_optimize(HeapSortAlg())
const TimSort   = maybe_optimize(TimSortAlg())
# Whenever InitialOptimizations is defined, RadixSort falls 
# back to Base.DEFAULT_STABLE which already incldues them.
const RadixSort = RadixSortAlg()

"""
    CombSort

Indicates that a sorting function should use the comb sort
algorithm. Comb sort traverses the collection multiple times
ordering pairs of elements with a given interval between them.
The interval decreases exponentially until it becomes 1, then
it switches to insertion sort on the whole input.

Characteristics:
 - *not stable* does not preserve the ordering of elements which
   compare equal (e.g. "a" and "A" in a sort of letters which
   ignores case).
 - *in-place* in memory.
 - *parallelizable* suitable for vectorization with SIMD instructions
   because it performs many independent comparisons.
 - *pathological inputs* such as `repeat(1:5.0, 4^8)` can make this algorithm perform very poorly.
 - *`n log n` average runtime* measured for random inputs of length up to 100 million, but theoretical runtime of `Θ(n^2)` for extremely long inputs.

## References
- Dobosiewicz, Wlodzimierz, (1980). "An efficient variation of bubble sort", Information Processing Letters, 11(1), pp. 5-6, https://doi.org/10.1016/0020-0190(80)90022-8.
 - Werneck, N. L., (2020). "ChipSort: a SIMD and cache-aware sorting module. JuliaCon Proceedings, 1(1), 12, https://doi.org/10.21105/jcon.00012
 - H. Inoue, T. Moriyama, H. Komatsu and T. Nakatani, "AA-Sort: A New Parallel Sorting Algorithm for Multi-Core SIMD Processors," 16th International Conference on Parallel Architecture and Compilation Techniques (PACT 2007), 2007, pp. 189-198, doi: 10.1109/PACT.2007.4336211.
"""
const CombSort  = maybe_optimize(CombSortAlg())

"""
    BranchyPatternDefeatingQuicksortAlg

Quicksort with improved performance on special input patterns.

Presorted inputs (including reverse and almost presorted ones), as well as inputs with many duplicates are
sorted in less than n log n time.
The code is based closely on the original C++ implementation by Orson Peters (see References).

Characteristics:
 - *not stable* does not preserve the ordering of elements which
   compare equal (e.g. "a" and "A" in a sort of letters which
   ignores case).
 - *in-place* in memory.
 - *`n log n` garuanteed runtime* by falling back to heapsort for pathological inputs.

## References
 - https://arxiv.org/pdf/2106.05123.pdf
 - https://github.com/orlp/pdqsort
"""
const BranchyPatternDefeatingQuicksort  = maybe_optimize(BranchyPatternDefeatingQuicksortAlg())
const BranchyPdqSort  = BranchyPatternDefeatingQuicksort

"""
    BranchlessPatternDefeatingQuicksortAlg

Quicksort with improved performance on special input patterns.

Presorted inputs (including reverse and almost presorted ones), as well as inputs with many duplicates are
sorted in less than n log n time. Uses branchless block partitioning scheme, which is faster for simple types.
The code is based closely on the original C++ implementation by Orson Peters (see References).

Characteristics:
 - *not stable* does not preserve the ordering of elements which
   compare equal (e.g. "a" and "A" in a sort of letters which
   ignores case).
 - *constant* auxilary memory (approximately 1KiB on 64-bit systems).
 - *`n log n` garuanteed runtime* by falling back to heapsort for pathological inputs.

## References
 - https://arxiv.org/pdf/2106.05123.pdf
 - https://github.com/orlp/pdqsort
 - https://dl.acm.org/doi/10.1145/3274660
 - http://arxiv.org/abs/1604.06697

"""
const BranchlessPatternDefeatingQuicksort  = maybe_optimize(BranchlessPatternDefeatingQuicksortAlg())
const BranchlessPdqSort  = BranchlessPatternDefeatingQuicksort

const PDQ_SMALL_THRESHOLD = 32
const PDQ_NINTHER_THRESHOLD = 128
const PDQ_PARTIAL_INSERTION_SORT_LIMIT = 8
const PDQ_BLOCK_SIZE = 64

## Heap sort

function sort!(v::AbstractVector, lo::Int, hi::Int, a::HeapSortAlg, o::Ordering)
    if lo > 1 || hi < length(v)
        return sort!(view(v, lo:hi), 1, length(v), a, o)
    end
    r = ReverseOrdering(o)
    heapify!(v, r)
    for i = length(v):-1:2
        # Swap the root with i, the last unsorted position
        x = v[i]
        v[i] = v[1]
        # The heap portion now ends at position i-1, but needs fixing up
        # starting with the root
        percolate_down!(v,1,x,r,i-1)
    end
    v
end


# Implementation of TimSort based on the algorithm description at:
#
#   http://svn.python.org/projects/python/trunk/Objects/listsort.txt
#   http://en.wikipedia.org/wiki/Timsort
#
# Original author: @kmsquire

const Run = UnitRange{Int}

const MIN_GALLOP = 7

mutable struct MergeState
    runs::Vector{Run}
    min_gallop::Int
end
MergeState() = MergeState(Run[], MIN_GALLOP)

# Determine a good minimum run size for efficient merging
# For details, see "Computing minrun" in
# http://svn.python.org/projects/python/trunk/Objects/listsort.txt
function merge_compute_minrun(N::Int, bits::Int)
    r = 0
    max_val = 2^bits
    while N >= max_val
        r |= (N & 1)
        N >>= 1
    end
    N + r
end
merge_compute_minrun(N::Int) = merge_compute_minrun(N, 6)

# Galloping binary search starting at left
# Finds the last value in v <= x
function gallop_last(o::Ordering, v::AbstractVector, x, lo::Int, hi::Int)
    i = lo
    inc = 1
    lo = lo-1
    hi = hi+1
    while i < hi && !lt(o, x, v[i])
        lo = i
        i += inc
        inc <<= 1
    end
    hi = min(i+1, hi)
    # Binary search
    while lo < hi-1
        i = (lo+hi)>>>1
        if lt(o, x, v[i])
            hi = i
        else
            lo = i
        end
    end
    lo
end

# Galloping binary search starting at right
# Finds the last value in v <= x
function rgallop_last(o::Ordering, v::AbstractVector, x, lo::Int, hi::Int)
    i = hi
    dec = 1
    lo = lo-1
    hi = hi+1
    while i > lo && lt(o, x, v[i])
        hi = i
        i -= dec
        dec <<= 1
    end
    lo = max(lo, i-1)
    # Binary search
    while lo < hi-1
        i = (lo+hi)>>>1
        if lt(o, x, v[i])
            hi = i
        else
            lo = i
        end
    end
    lo
end

# Galloping binary search starting at left
# Finds the first value in v >= x
function gallop_first(o::Ordering, v::AbstractVector, x, lo::Int, hi::Int)
    i = lo
    inc = 1
    lo = lo-1
    hi = hi+1
    while i < hi && lt(o, v[i], x)
        lo = i
        i += inc
        inc <<= 1
    end
    hi = min(i+1, hi)
    # Binary search
    while lo < hi-1
        i = (lo+hi)>>>1
        if lt(o, v[i], x)
            lo = i
        else
            hi = i
        end
    end
    hi
end

# Galloping binary search starting at right
# Finds the first value in v >= x
function rgallop_first(o::Ordering, v::AbstractVector, x, lo::Int, hi::Int)
    i = hi
    dec = 1
    lo = lo-1
    hi = hi+1
    while i > lo && !lt(o, v[i], x)
        hi = i
        i -= dec
        dec <<= 1
    end
    lo = max(lo, i-1)
    # Binary search
    while lo < hi-1
        i = (lo+hi)>>>1
        if lt(o, v[i], x)
            lo = i
        else
            hi = i
        end
    end
    hi
end

# Get the next run
# Returns the v range a:b, or b:-1:a for a reversed sequence
function next_run(o::Ordering, v::AbstractVector, lo::Int, hi::Int)
    lo == hi && return lo:hi
    if !lt(o, v[lo+1], v[lo])
        for i = lo+2:hi
            if lt(o, v[i], v[i-1])
                return lo:i-1
            end
        end
        return lo:hi
    else
        for i = lo+2:hi
            if !lt(o, v[i], v[i-1])
                return i-1:-1:lo
            end
        end
        return hi:-1:lo
    end
end

function merge_at(o::Ordering, v::AbstractVector, state::MergeState, n::Integer)
    a = state.runs[n]
    b = state.runs[n+1]
    merge(o,v,a,b,state)
    state.runs[n] = first(a):last(b)
    deleteat!(state.runs, n+1)
    nothing
end

# Merge consecutive runs
# For A,B,C,D = last four lengths, merge_collapse!()
# maintains 3 invariants:
#
#  A > B + C
#  B > C + D
#  C > D
#
# If any of these are violated, a merge occurs to
# correct it
function merge_collapse(o::Ordering, v::AbstractVector, state::MergeState)
    while true
        n = length(state.runs)
        n <= 1 && break

        # Check invariants 1 and 2
        if (n >= 3 && length(state.runs[end-2]) <= length(state.runs[end-1]) + length(state.runs[end])) ||
            (n >= 4 && length(state.runs[end-3]) <= length(state.runs[end-2]) + length(state.runs[end-1]))

            if length(state.runs[end-2]) < length(state.runs[end])
                merge_at(o,v,state,n-2)
            else
                merge_at(o,v,state,n-1)
            end

        # Check invariant 3
        elseif length(state.runs[end-1]) <= length(state.runs[end])
            merge_at(o,v,state,n-1)

        else # Invariant is satisfied
            break
        end
    end
end

# Merge runs a and b in vector v
function merge(o::Ordering, v::AbstractVector, a::Run, b::Run, state::MergeState)

    # First elements in a <= b[1] are already in place
    a = gallop_last(o, v, v[first(b)], first(a), last(a))+1: last(a)

    if length(a) == 0  return  end

    # Last elements in b >= a[end] are already in place
    b = first(b) : rgallop_first(o, v, v[last(a)], first(b), last(b))-1

    # Choose merge_lo or merge_hi based on the amount
    # of temporary memory needed (smaller of a and b)
    if length(a) < length(b)
        merge_lo(o, v, a, b, state)
    else
        merge_hi(o, v, a, b, state)
    end
end

# Merge runs a and b in vector v (a is smaller)
function merge_lo(o::Ordering, v::AbstractVector, a::Run, b::Run, state::MergeState)

    # Copy a
    v_a = v[a]

    # Pointer into (sub)arrays
    i = first(a)
    from_a = 1
    from_b = first(b)

    mode = :normal
    while true
        if mode == :normal
            # Compare and copy element by element
            count_a = count_b = 0
            while from_a <= length(a) && from_b <= last(b)
                if lt(o, v[from_b], v_a[from_a])
                    v[i] = v[from_b]
                    from_b += 1
                    count_a = 0
                    count_b += 1
                else
                    v[i] = v_a[from_a]
                    from_a += 1
                    count_a += 1
                    count_b = 0
                end
                i += 1

                # Switch to galloping mode if a string of elements
                # has come from the same set
                if count_b >= state.min_gallop || count_a >= state.min_gallop
                    mode = :galloping
                    break
                end
            end
            # Finalize if we've exited the loop normally
            if mode == :normal
                mode = :finalize
            end
        end

        if mode == :galloping
            # Use binary search to find range to copy
            while from_a <= length(a) && from_b <= last(b)
                # Copy the next run from b
                b_run = from_b : gallop_first(o, v, v_a[from_a], from_b, last(b)) - 1
                i_end = i + length(b_run) - 1
                v[i:i_end] = v[b_run]
                i = i_end + 1
                from_b = last(b_run) + 1

                # ... then copy the first element from a
                v[i] = v_a[from_a]
                i += 1
                from_a += 1

                if from_a > length(a) || from_b > last(b) break end

                # Copy the next run from a
                a_run = from_a : gallop_last(o, v_a, v[from_b], from_a, length(a))
                i_end = i + length(a_run) - 1
                v[i:i_end] = v_a[a_run]
                i = i_end + 1
                from_a = last(a_run) + 1

                # ... then copy the first element from b
                v[i] = v[from_b]
                i += 1
                from_b += 1

                # Return to normal mode if we haven't galloped...
                if length(a_run) < MIN_GALLOP && length(b_run) < MIN_GALLOP
                    mode = :normal
                    break
                end
                # Reduce min_gallop if this gallop was successful
                state.min_gallop -= 1
            end
            if mode == :galloping
                mode = :finalize
            end
            state.min_gallop = max(state.min_gallop,0) + 2  # penalty for leaving gallop mode
        end

        if mode == :finalize
            # copy end of a
            i_end = i + (length(a) - from_a)
            v[i:i_end] = v_a[from_a:end]
            break
        end
    end
end

# Merge runs a and b in vector v (b is smaller)
function merge_hi(o::Ordering, v::AbstractVector, a::Run, b::Run, state::MergeState)

    # Copy b
    v_b = v[b]

    # Pointer into (sub)arrays
    i = last(b)
    from_a = last(a)
    from_b = length(b)

    mode = :normal
    while true
        if mode == :normal
            # Compare and copy element by element
            count_a = count_b = 0
            while from_a >= first(a) && from_b >= 1
                if !lt(o, v_b[from_b], v[from_a])
                    v[i] = v_b[from_b]
                    from_b -= 1
                    count_a = 0
                    count_b += 1
                else
                    v[i] = v[from_a]
                    from_a -= 1
                    count_a += 1
                    count_b = 0
                end
                i -= 1

                # Switch to galloping mode if a string of elements
                # has come from the same set
                if count_b >= state.min_gallop || count_a >= state.min_gallop
                   mode = :galloping
                   break
                end
            end
            # Finalize if we've exited the loop normally
            if mode == :normal
                mode = :finalize
            end
        end

        if mode == :galloping
            # Use binary search to find range to copy
            while from_a >= first(a) && from_b >= 1
                # Copy the next run from b
                b_run = rgallop_first(o, v_b, v[from_a], 1, from_b) : from_b
                i_start = i - length(b_run) + 1
                v[i_start:i] = v_b[b_run]
                i = i_start - 1
                from_b = first(b_run) - 1

                # ... then copy the first element from a
                v[i] = v[from_a]
                i -= 1
                from_a -= 1

                if from_a < first(a) || from_b < 1 break end

                # Copy the next run from a
                a_run = rgallop_last(o, v, v_b[from_b], first(a), from_a) + 1: from_a
                i_start = i - length(a_run) + 1
                v[i_start:i] = v[a_run]
                i = i_start - 1
                from_a = first(a_run) - 1

                # ... then copy the first element from b
                v[i] = v_b[from_b]
                i -= 1
                from_b -= 1

                # Return to normal mode if we haven't galloped...
                if length(a_run) < MIN_GALLOP && length(b_run) < MIN_GALLOP
                    mode = :normal
                    break
                end
                # Reduce min_gallop if this gallop was successful
                state.min_gallop -= 1
            end
            if mode == :galloping
                mode = :finalize
            end
            state.min_gallop = max(state.min_gallop, 0) + 2  # penalty for leaving gallop mode
        end

        if mode == :finalize
            # copy start of b
            i_start = i - from_b + 1
            v[i_start:i] = v_b[1:from_b]
            break
        end
    end
end

# TimSort main method
function sort!(v::AbstractVector, lo::Int, hi::Int, ::TimSortAlg, o::Ordering)
    minrun = merge_compute_minrun(hi-lo+1)
    state = MergeState()
    i = lo
    while i <= hi
        run_range = next_run(o, v, i, hi)
        count = length(run_range)
        if count < minrun
            # Make a run of length minrun
            count = min(minrun, hi-i+1)
            run_range = i:i+count-1
            sort!(v, i, i+count-1, DEFAULT_STABLE, o)
        else
            if !issorted(run_range)
                run_range = last(run_range):first(run_range)
                reverse!(view(v, run_range))
            end
        end
        # Push this run onto the queue and merge if needed
        push!(state.runs, run_range)
        i = i+count
        merge_collapse(o, v, state)
    end
    # Force merge at the end
    while true
        n = length(state.runs)
        n <= 1 && break
        merge_at(o, v, state, n-1)
    end
    return v
end


function sort!(v::AbstractVector, lo::Int, hi::Int, ::CombSortAlg, o::Ordering)
    interval = (3 * (hi-lo+1)) >> 2

    while interval > 1
        @inbounds for j in lo:hi-interval
            a, b = v[j], v[j+interval]
            v[j], v[j+interval] = lt(o, b, a) ? (b, a) : (a, b)
        end
        interval = (3 * interval) >> 2
    end

    return sort!(v, lo, hi, InsertionSort, o)
end


## Radix sort
@static if VERSION >= v"1.9.0-DEV.482" # Base introduced radixsort in 1.9
    function sort!(vs::AbstractVector{T}, lo::Int, hi::Int, ::RadixSortAlg, o::Ordering, ts::Union{Nothing, AbstractVector{T}}=nothing) where T
        sort!(vs, lo, hi, Base.DEFAULT_STABLE, o)
    end
else

    # Map a bits-type to an unsigned int, maintaining sort order
    uint_mapping(::ForwardOrdering, x::Unsigned) = x
    for (signedty, unsignedty) in ((Int8, UInt8), (Int16, UInt16), (Int32, UInt32), (Int64, UInt64), (Int128, UInt128))
        # In Julia 0.4 we can just use unsigned() here
        @eval uint_mapping(::ForwardOrdering, x::$signedty) = reinterpret($unsignedty, xor(x, typemin(typeof(x))))
    end
    uint_mapping(::ForwardOrdering, x::Float32)  = (y = reinterpret(Int32, x); reinterpret(UInt32, ifelse(y < 0, ~y, xor(y, typemin(Int32)))))
    uint_mapping(::ForwardOrdering, x::Float64)  = (y = reinterpret(Int64, x); reinterpret(UInt64, ifelse(y < 0, ~y, xor(y, typemin(Int64)))))

    uint_mapping(::Sort.Float.Left, x::Float16)  = ~reinterpret(Int16, x)
    uint_mapping(::Sort.Float.Right, x::Float16)  = reinterpret(Int16, x)
    uint_mapping(::Sort.Float.Left, x::Float32)  = ~reinterpret(Int32, x)
    uint_mapping(::Sort.Float.Right, x::Float32)  = reinterpret(Int32, x)
    uint_mapping(::Sort.Float.Left, x::Float64)  = ~reinterpret(Int64, x)
    uint_mapping(::Sort.Float.Right, x::Float64)  = reinterpret(Int64, x)

    uint_mapping(rev::ReverseOrdering, x) = ~uint_mapping(rev.fwd, x)
    uint_mapping(::ReverseOrdering{ForwardOrdering}, x::Real) = ~uint_mapping(Forward, x) # maybe unnecessary; needs benchmark

    uint_mapping(o::By,   x     ) = uint_mapping(Forward, o.by(x))
    uint_mapping(o::Perm, i::Int) = uint_mapping(o.order, o.data[i])
    uint_mapping(o::Lt,   x     ) = error("uint_mapping does not work with general Lt Orderings")

    const RADIX_SIZE = 11
    const RADIX_MASK = 0x7FF

    function sort!(vs::AbstractVector{T}, lo::Int, hi::Int, ::RadixSortAlg, o::Ordering, ts::AbstractVector{T}=similar(vs)) where T
        # Input checking
        if lo >= hi;  return vs;  end

        # Make sure we're sorting a bits type
        OT = Base.Order.ordtype(o, vs)
        if !isbitstype(OT)
            error("Radix sort only sorts bits types (got $OT)")
        end

        # Init
        iters = ceil(Integer, sizeof(OT)*8/RADIX_SIZE)
        bin = zeros(UInt32, 2^RADIX_SIZE, iters)
        if lo > 1;  bin[1,:] .= lo-1;  end

        # Histogram for each element, radix
        for i = lo:hi
            v = uint_mapping(o, vs[i])
            for j = 1:iters
                idx = Int((v >> ((j-1)*RADIX_SIZE)) & RADIX_MASK) + 1
                @inbounds bin[idx,j] += 1
            end
        end

        # Sort!
        swaps = 0
        len = hi-lo+1
        for j = 1:iters
            # Unroll first data iteration, check for degenerate case
            v = uint_mapping(o, vs[hi])
            idx = Int((v >> ((j-1)*RADIX_SIZE)) & RADIX_MASK) + 1

            # are all values the same at this radix?
            if bin[idx,j] == len;  continue;  end

            cbin = cumsum(bin[:,j])
            ci = cbin[idx]
            ts[ci] = vs[hi]
            cbin[idx] -= 1

            # Finish the loop...
            @inbounds for i in hi-1:-1:lo
                v = uint_mapping(o, vs[i])
                idx = Int((v >> ((j-1)*RADIX_SIZE)) & RADIX_MASK) + 1
                ci = cbin[idx]
                ts[ci] = vs[i]
                cbin[idx] -= 1
            end
            vs,ts = ts,vs
            swaps += 1
        end

        if isodd(swaps)
            vs,ts = ts,vs
            for i = lo:hi
                vs[i] = ts[i]
            end
        end
        vs
    end

    function Base.Sort.Float.fpsort!(v::AbstractVector, ::RadixSortAlg, o::Ordering)
        @static if VERSION >= v"1.7.0-DEV"
            lo, hi = Base.Sort.Float.specials2end!(v, RadixSort, o)
        else
            lo, hi = Base.Sort.Float.nans2end!(v, o)
        end
        sort!(v, lo, hi, RadixSort, o)
    end
end

"""
    unguarded_insertion_sort!(v::AbstractVector, lo::Integer, hi::Integer, o::Ordering)

Sorts v[lo:hi] using insertion sort with the given ordering. Assumes
v[lo-1] is an element smaller than or equal to any element in v[lo:hi].
"""
function unguarded_insertion_sort!(v::AbstractVector, lo::Integer, hi::Integer, o::Ordering)
    lo_plus_1 = (lo + 1)::Integer
    @inbounds for i = lo_plus_1:hi
        j = i
        x = v[i]
        while true
            y = v[j-1]
            if !(lt(o, x, y)::Bool)
                break
            end
            v[j] = y
            j -= 1
        end
        v[j] = x
    end
    v
end

"""
    partial_insertion_sort!(v::AbstractVector, lo::Integer, hi::Integer, o::Ordering)

Attempts to use insertion sort on v[lo:hi]. Will return false if more than
PDQ_PARTIAL_INSERTION_SORT_LIMIT elements were moved, and abort sorting. Otherwise it will
successfully sort and return true.
"""
function partial_insertion_sort!(v::AbstractVector, lo::Integer, hi::Integer, o::Ordering)
    limit = 0
    lo_plus_1 = (lo + 1)::Integer
    @inbounds for i = lo_plus_1:hi
        j = i
        x = v[i]
        while j > lo
            y = v[j-1]
            if !(lt(o, x, y)::Bool)
                break
            end
            v[j] = y
            j -= 1
        end
        v[j] = x
        limit += i - j
        limit > PDQ_PARTIAL_INSERTION_SORT_LIMIT  && return false
    end
    return true
end

"""
    partition_right!(v::AbstractVector, lo::Integer, hi::Integer, a::BranchlessPatternDefeatingQuicksortAlg, o::Ordering, offsets_l::AbstractVector{Integer}, offsets_r::AbstractVector{Integer})

Partitions v[lo:hi] around pivot v[lo] using ordering o.

Elements equal to the pivot are put in the right-hand partition. Returns the position of the pivot
after partitioning and whether the passed sequence already was correctly partitioned. Assumes the
pivot is a median of at least 3 elements and that v[lo:hi] is at least PDQ_SMALL_THRESHOLD long.
Uses branchless partitioning.
"""
function partition_right!(v::AbstractVector, lo::Integer, hi::Integer, a::BranchlessPatternDefeatingQuicksortAlg, o::Ordering, offsets_l::AbstractVector{Int}, offsets_r::AbstractVector{Int})
    # input:
    # v[lo] -> pivot
    # output:
    # v[lo:pivot_position-1] < pivot
    # v[pivot_position] == pivot
    # v[pivot_position+1:hi] >= pivot
    @inbounds begin
        pivot = v[lo]

        # swap pointers
        # v[lo] is pivot -> start at lo + 1
        left = lo + 1
        right = hi
        # Find the first element greater than or equal than the pivot (the median of 3 guarantees
        # this exists).
        while lt(o, v[left], pivot)
            left += 1
        end
        # Find the first element strictly smaller than the pivot. We have to guard this search if
        # there was no element before v[left].
        if left - 1 == lo
            while left < right && !lt(o, v[right], pivot)
                right -= 1
            end
        else
            while !lt(o, v[right], pivot)
                right -= 1
            end
        end

        # If the first pair of elements that should be swapped to partition are the same element,
        # the passed in sequence already was correctly partitioned.
        was_already_partitioned = left >= right
        if !was_already_partitioned
            v[left], v[right] = v[right], v[left]
            left += 1
            right -= 1

            offsets_l_base = left
            offsets_r_base = right
            start_l = 0; start_r = 0
            num_l = 0; num_r = 0

            while left < right + 1
                # Fill up offset blocks with elements that are on the wrong side.
                # First we determine how much elements are considered for each offset block.
                num_unknown = right - left + 1
                left_split = num_l == 0 ? (num_r == 0 ? num_unknown ÷ 2 : num_unknown) : 0
                right_split = num_r == 0 ? (num_unknown - left_split) : 0

                # Fill the offset blocks.
                if left_split >= PDQ_BLOCK_SIZE
                    i = 0
                    while i < PDQ_BLOCK_SIZE
                        Cartesian.@nexprs 8 _ ->
                        begin
                            offsets_l[num_l+1] = i
                            num_l += Int(!lt(o, v[left], pivot))
                            left += 1
                            i += 1
                        end
                    end
                else
                    for i in 0:left_split-1
                        offsets_l[num_l+1] = i
                        num_l += Int(!lt(o, v[left], pivot))
                        left += 1
                    end
                end
                if right_split  >= PDQ_BLOCK_SIZE
                    i = 0
                    while i < PDQ_BLOCK_SIZE
                        Cartesian.@nexprs 8 _ ->
                        begin
                            offsets_r[num_r+1] = i
                            num_r += Int(lt(o, v[right], pivot))
                            right -= 1
                            i += 1
                        end
                    end
                else
                    for i in 0:right_split-1
                        offsets_r[num_r+1] = i
                        num_r += Int(lt(o, v[right], pivot))
                        right -= 1
                    end
                end

                # Swap elements and update block sizes and left/right boundaries.
                num = min(num_l, num_r)
                for i = 1:num
                    swap!(v, offsets_l_base + offsets_l[i+start_l], offsets_r_base - offsets_r[i+start_r])
                end
                num_l -= num; num_r -= num
                start_l += num; start_r += num

                if num_l == 0
                    start_l = 0
                    offsets_l_base = left
                end

                if num_r == 0
                    start_r = 0
                    offsets_r_base = right
                end
            end

            # We have now fully identified [left, right)'s proper position. Swap the last elements.
            if num_l > 0
                while num_l > 0
                    swap!(v, offsets_l_base + offsets_l[start_l+num_l], right)
                    num_l -= 1
                    right -= 1
                end
                left = right + 1
            end
            if num_r > 0
                while num_r > 0
                    swap!(v, left, offsets_r_base - offsets_r[start_r+num_r])
                    num_r -= 1
                    left += 1
                end
                right = left
            end

        end

        # Put the pivot in the right place.
        pivot_position = left - 1
        v[lo] = v[pivot_position]
        v[pivot_position] = pivot
    end
    return pivot_position, was_already_partitioned
end

"""
    partition_right!(v::AbstractVector, lo::Integer, hi::Integer, a::BranchyPatternDefeatingQuicksortAlg, o::Ordering, _, _)

Partitions v[lo:hi] around pivot v[lo] using ordering o.

Elements equal to the pivot are put in the right-hand partition. Returns the position of the pivot
after partitioning and whether the passed sequence already was correctly partitioned. Assumes the
pivot is a median of at least 3 elements and that v[lo:hi] is at least PDQ_SMALL_THRESHOLD long.
"""
function partition_right!(v::AbstractVector, lo::Integer, hi::Integer, a::BranchyPatternDefeatingQuicksortAlg, o::Ordering, _, _)
    # input:
    # v[lo] -> pivot
    # output:
    # v[lo:pivot_position-1] < pivot
    # v[pivot_position] == pivot
    # v[pivot_position+1:hi] >= pivot
    @inbounds begin
        pivot = v[lo]

        # swap pointers
        # v[lo] is pivot
        left = lo + 1
        right = hi
        # Find the left element greater than or equal than the pivot (the median of 3 guarantees
        # this exists).
        while lt(o, v[left], pivot)
            left += 1
        end
        # Find the first element strictly smaller than the pivot. We have to guard this search if
        # there was no element before v[left].
        if left - 1 == lo
            while left < right && !lt(o, v[right], pivot)
                right -= 1
            end
        else
            while !lt(o, v[right], pivot)
                right -= 1
            end
        end

        # If the first pair of elements that should be swapped to partition are the same element,
        # the passed in sequence already was correctly partitioned.
        was_already_partitioned = left >= right

        # Keep swapping pairs of elements that are on the wrong side of the pivot. Previously
        # swapped pairs guard the searches, which is why the first iteration is special-cased
        # above.
        while left < right
            swap!(v, left, right)
            left += 1
            right -= 1
            while lt(o, v[left], pivot)
                left += 1
            end
            while !lt(o, v[right], pivot)
                right -= 1
            end
        end

        # Put the pivot in the right place.
        pivot_position = left - 1
        v[lo] = v[pivot_position]
        v[pivot_position] = pivot

    end
    return pivot_position, was_already_partitioned
end

"""
    partition_left!(v::AbstractVector, lo::Integer, hi::Integer, o::Ordering)

Partitions v[lo:hi] around pivot v[lo] using ordering o.

Similar function to the one above, except elements equal to the pivot are put to the left of
the pivot and it doesn't check or return if the passed sequence already was partitioned.
Since this is rarely used (the many equal case), and in that case pdqsort already has O(n)
performance, no block quicksort is applied here for simplicity.
"""
function partition_left!(v::AbstractVector, lo::Integer, hi::Integer, o::Ordering)
    # input:
    # v[hi] -> pivot
    # output:
    # v[lo:pivot_position-1] <= pivot
    # v[pivot_position] == pivot
    # v[pivot_position+1:hi] > pivot
    
    @inbounds begin
        pivot = v[lo]
        left = lo + 1
        right = hi
        
        while lt(o, pivot, v[right])
            right -= 1
        end
        if right == hi
            while left < right && !lt(o, pivot, v[left])
                left += 1
            end
        else
            while !lt(o, pivot, v[left])
                left += 1
            end
        end
        
        while left < right
            swap!(v, left, right)
            while lt(o, pivot, v[right])
                right -= 1
            end
            while !lt(o, pivot, v[left])
                left += 1
            end
        end
        
        # Put the pivot in the right place.
        pivot_position = right
        v[lo] = v[pivot_position]
        v[pivot_position] = pivot
    end
    return pivot_position
end

# midpoint was added to Base.sort in version 1.4 and later moved to Base
# -> redefine for compatibility with earlier versions
_midpoint(lo::Integer, hi::Integer) = lo + ((hi - lo) >>> 0x01)

# modified from Base.sort
@inline function selectpivot!(v::AbstractVector, lo::Integer, hi::Integer, o::Ordering)
    @inbounds begin
        # use hi+1 to ensure reverse sorted list is swapped perfectly
        mi = _midpoint(lo, hi+1)
        # sort v[mi] <= v[lo] <= v[hi] such that the pivot is immediately in place
        if lt(o, v[lo], v[mi])
            v[mi], v[lo] = v[lo], v[mi]
        end
        
        if lt(o, v[hi], v[lo])
            if lt(o, v[hi],  v[mi])
                v[hi], v[lo], v[mi] = v[lo], v[mi], v[hi]
            else
                v[hi], v[lo] = v[lo], v[hi]
            end
        end
    end
end

@inline function swap!(v::AbstractVector, i::Integer, j::Integer)
    v[i], v[j] = v[j], v[i]
end

@inline function sort2!(v::AbstractVector, lo::Integer, hi::Integer, o::Ordering)
    lt(o, v[hi], v[lo]) && swap!(v, lo, hi)
end

@inline function sort3!(v::AbstractVector, lo::Integer, mid::Integer, hi::Integer, o::Ordering)
    sort2!(v,  lo, mid, o)
    sort2!(v, mid,  hi, o)
    sort2!(v,  lo, mid, o)
end

@inline function selectpivot_ninther(v::AbstractVector, lo::Integer, hi::Integer, len::Integer, o::Ordering)
    s2 = len ÷ 2
    sort3!(v, lo, lo + s2, hi, o)
    sort3!(v, lo + 1, lo + (s2 - 1), hi - 1, o)
    sort3!(v, lo + 2, lo + (s2 + 1), hi - 2, o)
    sort3!(v, lo + (s2 - 1), lo + s2, lo + (s2 + 1), o)
    swap!(v, lo, lo + s2)
end

pdqsort_loop!(v::AbstractVector, lo::Integer, hi::Integer, a::BranchlessPatternDefeatingQuicksortAlg, o::Ordering, bad_allowed::Integer, offsets_l::Nothing, offsets_r::Nothing, leftmost=true) =
pdqsort_loop!(v, lo, hi, a, o, bad_allowed, MVector{PDQ_BLOCK_SIZE, Int}(undef), MVector{PDQ_BLOCK_SIZE, Int}(undef), leftmost)

function pdqsort_loop!(v::AbstractVector, lo::Integer, hi::Integer, a::PatternDefeatingQuicksortAlg, o::Ordering, bad_allowed::Integer, offsets_l, offsets_r, leftmost=true)
    # Use a while loop for tail recursion elimination.
    @inbounds while true
        len = hi - lo + 1
        # Insertion sort is faster for small arrays.
        if len < PDQ_SMALL_THRESHOLD
            if leftmost
                sort!(v, lo, hi, InsertionSort, o)
            else
                unguarded_insertion_sort!(v, lo, hi, o)
            end
            return v
        end
        
        # Choose pivot as median of 3 or pseudomedian of 9.
        if len > PDQ_NINTHER_THRESHOLD
            selectpivot_ninther(v, lo, hi, len, o)
        else
            selectpivot!(v, lo, hi, o)
        end
        # If v[lo - 1] is the end of the right partition of a previous partition operation
        # there is no element in [begin, end) that is smaller than v[lo - 1]. Then if our
        # pivot compares equal to v[lo - 1] we change strategy, putting equal elements in
        # the left partition, greater elements in the right partition. We do not have to
        # recurse on the left partition, since it's sorted (all equal).
        if !leftmost && !lt(o, v[lo-1], v[lo])
            lo = partition_left!(v, lo, hi, o) + 1
            continue
        end
        
        # Partition and get results.
        pivot_pos, was_already_partitioned = partition_right!(v, lo, hi, a, o, offsets_l, offsets_r)
        
        # Check for a highly unbalanced partition.
        l_len = pivot_pos - lo;
        r_len = hi - (pivot_pos + 1);
        is_highly_unbalanced = l_len < len ÷ 8 || r_len < len ÷ 8
        
        # If we got a highly unbalanced partition we shuffle elements to break many patterns.
        if is_highly_unbalanced
            # If we had too many bad partitions, switch to heapsort to guarantee O(n log n).
            bad_allowed -= 1
            if bad_allowed <= 0
                sort!(v, lo, hi, HeapSortAlg(), o)
                return v
            end
            
            if l_len > PDQ_SMALL_THRESHOLD
                swap!(v,    lo,             lo + l_len ÷ 4)
                swap!(v,    pivot_pos - 1,  pivot_pos - l_len ÷ 4)
                
                if (l_len > PDQ_NINTHER_THRESHOLD)
                    swap!(v, lo + 1,        lo + (l_len ÷ 4 + 1))
                    swap!(v, lo + 2,        lo + (l_len ÷ 4 + 2))
                    swap!(v, pivot_pos - 2, pivot_pos - (l_len ÷ 4 + 1))
                    swap!(v, pivot_pos - 3, pivot_pos - (l_len ÷ 4 + 2))
                end
            end

            if r_len > PDQ_SMALL_THRESHOLD
                swap!(v,    pivot_pos + 1,  pivot_pos + (1 + r_len ÷ 4))
                swap!(v,    hi,             hi - r_len ÷ 4)
                
                if (r_len > PDQ_NINTHER_THRESHOLD)
                    swap!(v, pivot_pos + 2, pivot_pos + (2 + r_len ÷ 4))
                    swap!(v, pivot_pos + 3, pivot_pos + (3 + r_len ÷ 4))
                    swap!(v, hi - 1 ,       hi - 1 - r_len ÷ 4)
                    swap!(v, hi - 2,        hi - 2 - r_len ÷ 4)
                end
            end
        else
            # If we were decently balanced and we tried to sort an already partitioned
            # sequence try to use insertion sort.
            if was_already_partitioned &&
                partial_insertion_sort!(v, lo, pivot_pos, o) &&
                partial_insertion_sort!(v, pivot_pos + 1, hi, o)
                return v
            end
        end
        
        # Sort the left partition first using recursion and do tail recursion elimination for
        # the right-hand partition.
        pdqsort_loop!(v, lo, pivot_pos-1, a, o, bad_allowed, offsets_l, offsets_r, leftmost)
        lo = pivot_pos + 1
        leftmost = false
    end
end

# integer logarithm base two, ignoring sign
function log2i(n::Integer)
    sizeof(n) << 3 - leading_zeros(abs(n))
end

sort!(v::AbstractVector, lo::Int, hi::Int, a::PatternDefeatingQuicksortAlg, o::Ordering) =
pdqsort_loop!(v, lo, hi, a, o, log2i(hi + 1 - lo), nothing, nothing)

#=
This implementation of pattern-defeating quicksort is based on the original code from Orson Peters,
available at https://github.com/orlp/pdqsort.
Original license notice:
"""
Copyright (c) 2021 Orson Peters <orsonpeters@gmail.com>

This software is provided 'as-is', without any express or implied warranty. In no event will the
authors be held liable for any damages arising from the use of this software.

Permission is granted to anyone to use this software for any purpose, including commercial
applications, and to alter it and redistribute it freely, subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the
   original software. If you use this software in a product, an acknowledgment in the product
   documentation would be appreciated but is not required.

2. Altered source versions must be plainly marked as such, and must not be misrepresented as
   being the original software.

3. This notice may not be removed or altered from any source distribution.
"""
=#
end # module
