"Value and index tuple: enables sortandperm_radix"
Valindex{T, S<:Integer} = Tuple{T, S}

isbits(::Type{Valindex{T,S}}) where {T,S} = Base.isbits(T)
uint_mapping(o, vi::Valindex{T,S}) where {T,S} = uint_mapping(o, vi[1]) # enable sorting
ordtype(o, vs::Valindex{T,S}) where {T,S} = Base.ordtype(o, vs[1])
sizeof(::Type{Valindex{T,S}}) where {T,S} = Base.sizeof(T)

"""
	sortandperm(v, alg, o)

returns both the sort(v) as well as sortperm(v)
"""
function sortandperm(v, alg::RadixSortAlg, o::Ordering= Base.ForwardOrdering())
	sortandperm_radix(v, o)
end

function _sortandperm_radix(v::AbstractVector{T}, o::Ordering= Base.ForwardOrdering()) where T
	vv = Valindex{T,Int}[(vv,i) for (i,vv) in enumerate(v)]
	sort!(vv, alg=RadixSort, order = o)
	vv
end

function sortandperm_radix(v::AbstractVector{T}, o::Ordering= Base.ForwardOrdering()) where T
	vv = _sortandperm_radix(v, o)
	res = Int[v[2] for v in vv]
	val = T[v[1] for v in vv]
	return (val, res)
end

"""
	sortperm_radix(v, o)

sortperm using the radixsort algorithm
"""
function sortperm_radix(v, o::Ordering = Base.ForwardOrdering())
	vv = _sortandperm_radix(v,o)
	Int[v[2] for v in vv]
end