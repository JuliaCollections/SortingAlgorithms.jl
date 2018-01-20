__precompile__()

module SortingAlgorithms

using Compat
using DataStructures
import Compat.view
using Base.Sort
using Base.Order

import Base.Sort: sort!
import DataStructures: heapify!, percolate_down!

export HeapSort, TimSort, RadixSort

struct HeapSortAlg  <: Algorithm end
struct TimSortAlg   <: Algorithm end
struct RadixSortAlg <: Algorithm end

const HeapSort  = HeapSortAlg()
const TimSort   = TimSortAlg()
const RadixSort = RadixSortAlg()


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


## Radix sort

# Map a bits-type to an unsigned int, maintaining sort order
uint_mapping(::ForwardOrdering, x::Unsigned) = x
for (signedty, unsignedty) in ((Int8, UInt8), (Int16, UInt16), (Int32, UInt32), (Int64, UInt64), (Int128, UInt128))
    # In Julia 0.4 we can just use unsigned() here
    @eval uint_mapping(::ForwardOrdering, x::$signedty) = reinterpret($unsignedty, xor(x, typemin(typeof(x))))
end
uint_mapping(::ForwardOrdering, x::Float32)  = (y = reinterpret(Int32, x); reinterpret(UInt32, ifelse(y < 0, ~y, xor(y, typemin(Int32)))))
uint_mapping(::ForwardOrdering, x::Float64)  = (y = reinterpret(Int64, x); reinterpret(UInt64, ifelse(y < 0, ~y, xor(y, typemin(Int64)))))

uint_mapping(rev::ReverseOrdering, x) = ~uint_mapping(rev.fwd, x)
uint_mapping(::ReverseOrdering{ForwardOrdering}, x::Real) = ~uint_mapping(Forward, x) # maybe unnecessary; needs benchmark

uint_mapping(o::By,   x     ) = uint_mapping(Forward, o.by(x))
uint_mapping(o::Perm, i::Int) = uint_mapping(o.order, o.data[i])
uint_mapping(o::Lt,   x     ) = error("uint_mapping does not work with general Lt Orderings")

const RADIX_SIZE = 11
const RADIX_MASK = 0x7FF

function sort!(vs::AbstractVector, lo::Int, hi::Int, ::RadixSortAlg, o::Ordering, ts=similar(vs))
    # Input checking
    if lo >= hi;  return vs;  end

    # Make sure we're sorting a bits type
    T = Base.Order.ordtype(o, vs)
    if !isbits(T)
        error("Radix sort only sorts bits types (got $T)")
    end

    # Init
    iters = ceil(Integer, sizeof(T)*8/RADIX_SIZE)
    bin = zeros(UInt32, 2^RADIX_SIZE, iters)
    if lo > 1;  bin[1,:] = lo-1;  end

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
    lo, hi = Base.Sort.Float.nans2end!(v,o)
    sort!(v, lo, hi, RadixSort, o)
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

end # module
