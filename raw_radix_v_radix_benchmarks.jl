using SortMark, SortingAlgorithms
using Base.Sort
using Base.Order

struct RawRadixSort2Alg <: Algorithm end
const RawRadixSort2 = RawRadixSort2Alg()
function Base.sort!(v::AbstractVector{<:Unsigned}, lo::Integer, hi::Integer, ::RawRadixSort2Alg, o::ForwardOrdering)
    #compression, bits, chunk_size = heuristic(typemin(eltype(v)), typemax(eltype(v)), hi-lo+1)
    out = SortingAlgorithms.radix_sort!(v, similar(v), lo, hi, unsigned(8*sizeof(eltype(v))), Val(0xb))
    out === v || copyto!(v, out)
    v
end

for T in sort(SortMark.UInts, by=sizeof)
    df = SortMark.make_df([RawRadixSort2, RadixSort],
        Types=[T], orders=[Base.Order.Forward],
        lens=SortMark.lengths(2, 2_000_000, 5),
        sources=Dict(:simple=>SortMark.sources[:simple]),
        seconds=nothing, samples=20)
    compute!(df)
    stat!(df)
    display(df[:, [:len, :Type, :pvalue, :point_estimate, :confint]])

    plot(log10.(df.len), first.(df.confint),
        title="$T: 95% CI of runtime ratio: This/That",
        xlabel="log10(length)",
        legend=false)
    display(plot!(log10.(df.len), last.(df.confint)))
end
