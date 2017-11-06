import SortingAlgorithms: uint_mapping, RADIX_SIZE, RADIX_MASK, RadixSortAlg
using Compat


import Base: sortperm!, Ordering, Order

lo = 1
hi = length(vs)
o = Base.ForwardOrdering()


function sortperm!(vs::AbstractVector, lo::Int, hi::Int, ::RadixSortAlg, o::Ordering, ts=similar(vs))
    # Input checking
	print("ehllo")
	if lo >= hi;  return vs;  end

	orderindex = zeros(Int32, length(vs))
	orderindex1 = zeros(Int32, length(vs))

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
	        idx = @compat(Int((v >> (j-1)*RADIX_SIZE) & RADIX_MASK)) + 1
	        @inbounds bin[idx,j] += 1
	    end
	end

	# Sort!
	swaps = 0
	len = hi-lo+1
	for j = 1:iters
	    # Unroll first data iteration, check for degenerate case
	    v = uint_mapping(o, vs[hi])
	    idx = @compat(Int((v >> (j-1)*RADIX_SIZE) & RADIX_MASK)) + 1

	    # are all values the same at this radix?

	    if bin[idx,j] == len;  continue;  end
	    cbin = cumsum(bin[:,j])
	    ci = cbin[idx]
	    ts[ci] = vs[hi]
	    orderindex[ci] = hi
	    cbin[idx] -= 1

	    # Finish the loop...
	    @inbounds for i in hi-1:-1:lo
	        v = uint_mapping(o, vs[i])
	        idx = @compat(Int((v >> (j-1)*RADIX_SIZE) & RADIX_MASK)) + 1
	        ci = cbin[idx]
	        ts[ci] = vs[i]
	        orderindex[ci] = i
	        cbin[idx] -= 1
	    end
	    vs,ts = ts,vs
	    orderindex, orderindex1 = orderindex1, orderindex
	    swaps += 1
	end

	#return (orderindex,orderindex1)

	if isodd(swaps)
		return orderindex1
	else
		return orderindex
	end
end

using StatsBase
vs = sample(Int32(1):Int32(10), 10, replace = false)
ts=similar(vs)

vs
cvs = copy(vs)
abc = sortperm!(cvs,lo, hi,RadixSortAlg(),o)

vs[abc]
