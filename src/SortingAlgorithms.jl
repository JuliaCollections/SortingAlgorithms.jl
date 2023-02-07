__precompile__()

module SortingAlgorithms

using DataStructures
using Base.Sort
using Base.Order

import Base.Sort: sort!
import DataStructures: heapify!, percolate_down!
import StaticArrays: MVector

export HeapSort, TimSort, RadixSort, CombSort, PagedMergeSort, ThreadedPagedMergeSort

struct HeapSortAlg  <: Algorithm end
struct TimSortAlg   <: Algorithm end
struct RadixSortAlg <: Algorithm end
struct CombSortAlg  <: Algorithm end
struct PagedMergeSortAlg  <: Algorithm end
struct ThreadedPagedMergeSortAlg  <: Algorithm end

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
    PagedMergeSort

Indicates that a sorting function should use the paged merge sort
algorithm. Paged merge sort uses is a merge sort, that uses different
merge routines to achieve stable sorting with a scratch space of size O(√n).
The merge routine for merging large subarrays merges
blocks/pages of size O(√n) almost in place, before reordering them using a page table.
At deeper recursion levels, where the scratch space is big enough,
normal merging is used, where one input is copied into the scratch space.
When the scratch space is large enough to hold the complete subarray,
the input is merged interleaved from both sides, which increases performance
for random data.

Characteristics:
 - *stable* does preserve the ordering of elements which
   compare equal (e.g. "a" and "A" in a sort of letters which
   ignores case).
 - *O(√n)* auxilary memory usage.
 - *`O(n log n)` garuanteed runtime*.

## References
 - https://link.springer.com/chapter/10.1007/BFb0016253
 - https://max-arbuzov.blogspot.com/2021/10/merge-sort-with-osqrtn-auxiliary-memory.html
"""
const PagedMergeSort  = PagedMergeSortAlg()

"""
    ThreadedPagedMergeSort

Multithreaded version of PagedMergeSort using Threads.nthreads-times the auxilary space.
Uses multithreaded recursion (not multithreaded merging), so the maximum speedup is
limited to O(log n)
Characteristics:
 - *stable* does preserve the ordering of elements which
   compare equal (e.g. "a" and "A" in a sort of letters which
   ignores case).
 - *O(√n)* auxilary memory usage.
 - *`O(n log n)` garuanteed runtime*.

## References
 - https://link.springer.com/chapter/10.1007/BFb0016253
 - https://max-arbuzov.blogspot.com/2021/10/merge-sort-with-osqrtn-auxiliary-memory.html
 - https://en.wikipedia.org/wiki/Merge_sort#Merge_sort_with_parallel_recursion
"""
const ThreadedPagedMergeSort  = ThreadedPagedMergeSortAlg()

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

###
# ThreadedPagedMergeSort
###

# merge v[lo:hiA] and v[hiA+1:hi] ([A;B])  using buffer t[1:1 + hi-lo]
function twoended_merge!(v::AbstractVector{T}, t::AbstractVector{T}, lo::Integer, hi::Integer, hiA::Integer, o::Ordering) where T
    @assert lo <= hiA <= hi
    loA = lo
    loB = hiA + 1
    hiB = hi
    
    # output array indices 
    oL = 1
    oR = 1 + hi-lo
    
    # input array indices
    iAL = loA
    iBL = loB
    iAR = hiA
    iBR = hiB
    
    @inbounds begin
        # two ended merge
        while iAL < iAR && iBL < iBR
            if lt(o,v[iBL], v[iAL])
                t[oL] = v[iBL]
                iBL += 1
            else
                t[oL] = v[iAL]
                iAL += 1
            end
            oL +=1
            
            if lt(o,v[iAR], v[iBR])
                t[oR] = v[iBR]
                iBR -= 1
            else
                t[oR] = v[iAR]
                iAR -= 1
            end
            oR -=1
        end        
        # cleanup
        # regular merge
        while iAL <= iAR && iBL <= iBR
            if  lt(o,v[iBL], v[iAL])
                t[oL] = v[iBL]
                iBL += 1
            else
                t[oL] = v[iAL]
                iAL += 1
            end
            oL += 1
        end
        # either A or B is empty -> copy remaining items
        while iAL <= iAR
            t[oL] = v[iAL]
            iAL += 1
            oL += 1
        end
        while iBL <= iBR
            t[oL] = v[iBL]
            iBL += 1
            oL += 1
        end
        # copy back from t to v
        offset = lo-1
        len = 1 + hi - lo
        @inbounds for i = 1:len
            v[offset+i] = t[i]
        end
    end
end

# merge v[lo:lo+lenA-1] and v[lo+lenA:hi] using buffer t[1:lenA]
# based on Base.Sort MergeSort
function merge!(v::AbstractVector{T},t::AbstractVector{T}, lenA::Integer, lo::Integer, hi::Integer, o::Ordering) where T
    @inbounds begin
        i = 1
        j = lo
        while i <= lenA
            t[i] = v[j]
            i +=1
            j +=1
        end
        iA = 1
        k = lo
        iB = lo + lenA
        while k < iB <= hi
            if lt(o,v[iB], t[iA])
                v[k] = v[iB]
                iB += 1
            else
                v[k] = t[iA]
                iA += 1
            end
            k += 1
        end
        while iA <= lenA
            v[k] = t[iA]
            k += 1
            iA += 1
        end
    end
end

# macro used for block management in pagedMerge!
macro getNextBlock!()
    quote
        if iA > nextBlockA * blocksize + lo
            currentBlock = nextBlockA
            nextBlockA += 1                
        else
            currentBlock = nextBlockB
            nextBlockB += 1
        end
        blockLocation[currentBlockIdx] = currentBlock
        currentBlockIdx += 1
    end |> esc
end

# merge v[lo:endA] and v[endA+1:hi] using buffer buf in O(sqrt(n)) space
function pagedMerge!(v::AbstractVector{T}, buf::AbstractVector{T}, lo::Integer, endA::Integer, hi::Integer, blockLocation::AbstractVector{<:Integer}, o::Ordering) where T
    @assert lo < endA < hi
    iA = lo
    iB = endA + 1
    endB = hi    
    lenA = endA + 1 - lo
    lenB = endB - endA

    # regular merge if buffer is big enough
    if lenA <= length(buf)
        merge!(v,buf,lenA,lo,hi,o)
        return
    elseif lenB <= length(buf)
        # TODO ?
        # does not occur in balanced mergesort where length(A) <= length(B)
        error("not implemented")
        return
    end

    len = hi + 1 - lo
    blocksize = isqrt(len)
    nBlocks = len ÷ blocksize
    @assert length(buf) >= 3blocksize
    @assert length(blockLocation) >= nBlocks+1

    @inline getBlockOffset(block) = (block-1)*blocksize + lo - 1         

    @inbounds begin 
        ##################
        # merge
        ##################
        # merge into buf until full
        oBuf = 1
        while oBuf <= 3blocksize    # cannot run out of input elements here
            if lt(o, v[iB], v[iA])  # -> merge! would have been used
                buf[oBuf] = v[iB]
                iB += 1
            else
                buf[oBuf] = v[iA]
                iA += 1
            end
            oBuf += 1
        end

        nextBlockA = 1
        nextBlockB = (endA+blocksize-lo) ÷ blocksize + 1
        blockLocation .= 0
        blockLocation[1:3] = -1:-1:-3

        oIdx = 1
        currentBlock = 0
        currentBlockIdx = 4
        # more efficient loop while more than blocksize elements of A and B are remaining
        while iA < endA-blocksize && iB < endB-blocksize
            @getNextBlock!
            offset = (currentBlock-1)*blocksize
            oIdx = lo + offset
            while oIdx <= blocksize+offset + lo - 1
                if lt(o, v[iB], v[iA])
                    v[oIdx] = v[iB]
                    iB += 1
                else
                    v[oIdx] = v[iA]
                    iA += 1
                end
                oIdx += 1
            end
        end
        # merge until either A or B is empty
        while iA <= endA && iB <= endB
            @getNextBlock!
            oIdx = 1
            offset = getBlockOffset(currentBlock)
            while oIdx <= blocksize && iA <= endA && iB <= endB
                if lt(o, v[iB], v[iA])
                    v[offset+oIdx] = v[iB]
                    iB += 1
                else
                    v[offset+oIdx] = v[iA]
                    iA += 1
                end
                oIdx += 1
            end
        end
        # copy remaining elements
        # either A or B is empty
        # copy rest of A
        while iA <= endA
            if oIdx > blocksize
                @getNextBlock!
                oIdx = 1
            end
            offset = getBlockOffset(currentBlock)
            while oIdx <= blocksize && iA <= endA
                v[offset + oIdx] = v[iA]
                iA += 1
                oIdx += 1
            end
        end
        # copy rest of B
        while iB <= endB
            if oIdx > blocksize
                @getNextBlock!
                oIdx = 1
            end
            offset = getBlockOffset(currentBlock)
            while oIdx <= blocksize && iB <= endB
                v[offset + oIdx] = v[iB]
                iB += 1
                oIdx += 1
            end
        end
        # copy last partial block to end
        partialBlockPresent = oIdx <= blocksize
        if partialBlockPresent
            offset = getBlockOffset(currentBlock)
            offset2 = nBlocks*blocksize + lo - 1
            for j = 1:oIdx-1
                v[offset2 + j] = v[offset + j]
            end
            blockLocation[currentBlockIdx-1] = 0
        end
        #########################################
        # calculate location of the 3 free blocks
        #########################################
        nFreeBlocksB = nBlocks + 1 - nextBlockB
        nFreeBlocksA = 3 - nFreeBlocksB - Int(partialBlockPresent)
        freeBlocks = MVector{3,Int}(undef)
        i = 1
        for j = 0:nFreeBlocksA-1
            freeBlocks[i] = nextBlockA + j
            i += 1
        end
        for j = 0:nFreeBlocksB-1
            freeBlocks[i] = nextBlockB + j
            i += 1
        end
        if partialBlockPresent
            freeBlocks[i] = currentBlock
        end       
        freeBlocksIdx = 3
        doneBlockIdx = 1
        currentBlock = freeBlocks[end]
        ##################
        # rearrange blocks
        ##################
        while true
            blc = blockLocation[currentBlock] # index of block with data belonging to currentBlock
            if blc > 0
                # data for currentBlock is in v                
                offset = getBlockOffset(currentBlock)
                offset2 = getBlockOffset(blc)
                for j = 1:blocksize
                    v[offset + j] = v[offset2 + j]
                end
                blockLocation[currentBlock] = 0
                currentBlock = blc
            else
                # data for currentBlock is in buf
                offset = getBlockOffset(currentBlock)
                offset2 = (-blc-1)*blocksize
                for j = 1:blocksize
                    v[offset + j] = buf[offset2 + j]
                end
                blockLocation[currentBlock] = 0
                if freeBlocksIdx > 1
                    # get next free block
                    freeBlocksIdx -= 1
                    currentBlock = freeBlocks[freeBlocksIdx]
                else
                    # no free block remains
                    # make sure that all blocks are done
                    while blockLocation[doneBlockIdx] == 0 || blockLocation[doneBlockIdx] == doneBlockIdx
                        doneBlockIdx += 1
                        doneBlockIdx == nBlocks && return
                    end
                    # copy misplaced block into buf and continue        
                    currentBlock = blockLocation[doneBlockIdx]
                    offset = getBlockOffset(currentBlock)
                    for j = 1:blocksize
                        buf[j] = v[offset + j]
                    end
                    blockLocation[doneBlockIdx] = -1
                end
            end
        end
    end    
end

function pagedmergesort!(v::AbstractVector{T}, lo::Integer, hi::Integer, buf::AbstractVector{T}, blockLocation, o=Base.Order.Forward) where T
    len = hi + 1 -lo
    if len <= Base.SMALL_THRESHOLD
        return Base.Sort.sort!(v, lo, hi, Base.Sort.InsertionSortAlg(), o)
    end
    m = Base.midpoint(lo,hi)
    pagedmergesort!(v,lo,m,buf,blockLocation,o)
    pagedmergesort!(v,m+1,hi,buf,blockLocation,o)
    if len <= length(buf)
        twoended_merge!(v, buf, lo, hi, m,o)
    else
        pagedMerge!(v, buf, lo, m, hi, blockLocation, o)        
    end
    return v
end

function threaded_pagedmergesort!(v::AbstractVector, lo::Integer, hi::Integer, bufs, blockLocations, c::Channel, threadingThreshold::Integer, o=Base.Order.Forward)       
    len = hi + 1 -lo
    if len <= Base.SMALL_THRESHOLD
        return Base.Sort.sort!(v, lo, hi, Base.Sort.InsertionSortAlg(), o)
    end
    m = Base.midpoint(lo,hi)
    if len > threadingThreshold
        thr = Threads.@spawn threaded_pagedmergesort!(v,lo,m,bufs,blockLocations,c,threadingThreshold,o)
        threaded_pagedmergesort!(v,m+1,hi,bufs,blockLocations,c,threadingThreshold,o)
        wait(thr)
        id = take!(c)
        buf = bufs[id]
        blockLocation = blockLocations[id]
    else
        id = take!(c)
        buf = bufs[id]
        blockLocation = blockLocations[id]
        pagedmergesort!(v,lo,m,buf,blockLocation,o)
        pagedmergesort!(v,m+1,hi,buf,blockLocation,o)
    end
    if len <= length(buf)
        twoended_merge!(v, buf, lo, hi, m, o)
    else
        pagedMerge!(v, buf, lo, m, hi, blockLocation, o)        
    end
    put!(c,id)
    return v
end

const PAGEDMERGESORT_THREADING_THRESHOLD = 2^13

function sort!(v::AbstractVector, lo::Integer, hi::Integer, a::PagedMergeSortAlg, o::Ordering)
    lo >= hi && return v
    n = hi + 1 - lo
    blocksize = isqrt(n)
    buf = Vector{eltype(v)}(undef,3blocksize)
    nBlocks = n ÷ blocksize
    blockLocation = Vector{Int}(undef,nBlocks+1)
    pagedmergesort!(v,lo,hi,buf,blockLocation,o)
    return v
end

function sort!(v::AbstractVector, lo::Integer, hi::Integer, a::ThreadedPagedMergeSortAlg, o::Ordering)
    lo >= hi && return v
    n = hi + 1 - lo
    nThreads=Threads.nthreads()
    (n < PAGEDMERGESORT_THREADING_THRESHOLD || nThreads < 2) && return sort!(v, lo, hi, PagedMergeSort, o)
    threadingThreshold = max(n ÷ 4nThreads, PAGEDMERGESORT_THREADING_THRESHOLD)
    blocksize = isqrt(n)
    nBlocks = n ÷ blocksize
    bufs = [Vector{eltype(v)}(undef,3blocksize) for _ in 1:nThreads] # allocate buffer for each thread
    blockLocation = [Vector{Int}(undef,nBlocks+1) for _ in 1:nThreads]
    c = Channel{Int}(nThreads) # channel holds indices of available buffers
    for i=1:nThreads
        put!(c,i)
    end
    threaded_pagedmergesort!(v,lo,hi,bufs,blockLocation,c,threadingThreshold,o)
    return v
end
end # module
