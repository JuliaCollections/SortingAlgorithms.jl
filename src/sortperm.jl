using SortingAlgorithms
import SortingAlgorithms: uint_mapping, RADIX_SIZE, RADIX_MASK, RadixSortAlg
using Compat, Base.Test
import Base: sortperm, sortperm!, Ordering, Order, setindex!, isbits, ordtype, sizeof

"Value and index tuple: enables sortandperm_radix"
Valindex{T, S<:Integer} = Tuple{T, S}

isbits(::Type{Valindex{T,S}}) where {T,S} = isbits(T)
uint_mapping(o, vi::Valindex{T,S}) where {T,S} = uint_mapping(o, vi[1]) # enable sorting
ordtype(o, vs::Valindex{T,S}) where {T,S} = ordtype(o, vs[1])
sizeof(::Type{Valindex{T,S}}) where {T,S} = sizeof(T)

"returns both the sort(v) as well as sortperm(v)"
function sortandperm(v, alg::RadixSortAlg)
	sortandperm_radix(v)
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

function sortperm_radix(v, o::Ordering = Base.ForwardOrdering())
	vv = _sortandperm_radix(v,o)
	Int[v[2] for v in vv]
end

"""RadixSort perm"""
function sortperm(v::AbstractVector,
                  alg::RadixSortAlg;
                  lt=isless,
                  by=identity,
                  rev::Bool = false,
                  order::Ordering=Base.Forward)
    ordr = Base.ord(lt,by,rev,order)
    vv = _sortandperm_radix(v, ordr)
	Int[v[2] for v in vv]
end

function sortperm!(x::AbstractVector{<:Integer}, v::AbstractVector,
                   alg::RadixSortAlg;
                   lt=isless,
                   by=identity,
                   rev::Bool = false,
                   order::Ordering=Base.Forward,
                   initialized::Bool=false)
    ordr = Base.ord(lt,by,rev,order)

    vv = _sortandperm_radix(v, ordr)

	for i = 1:length(vv)
		x[i] = vv[i][2]
	end
	x
end
