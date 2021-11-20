using SortMark
using SortingAlgorithms
using SortingAlgorithms: Serializable, serialize, heuristic, compress!, radix_sort!, deserialize!, time_est
using Base.Sort
import Base.Sort.sort!
using Base.Order
using Optim
using Random: shuffle!
using Statistics
using Plots

## Esitmate
data = []
    #=Juno.@profiler=# for T in sort(SortMark.UInts, by=sizeof)
    for bits in 1:(8*sizeof(T))
        print("\r",bits,"/$(8*sizeof(T))")
        lens = SortMark.lengths(30,3_000_000÷sizeof(T),2)
        for len in lens
            if rand() < .1
                x = rand(T(0):(T(1)<<bits)-T(1), len)
                max_chunk = ceil(log2(len))+1
                mins = []
                for chunk in 1:max_chunk
                    times = []
                    end_time = time() + 1/(8*sizeof(T))/length(lens)/max_chunk
                    while true
                        c = copy(x)
                        ts = similar(c)
                        t = @elapsed SortingAlgorithms.radix_sort!(c, ts,
                            firstindex(c), lastindex(c), unsigned(bits), Val(UInt8(chunk)))
                        push!(times, t)
                        length(times) < 10 || time() < end_time || break
                    end
                    push!(mins, minimum(times))
                end
                push!(data, ((T, bits, len), mins, minimum(mins)))
            end
        end
    end
end

function heuristic_k(T, bits, length, k)
    1#max(1, round(Int, log(length)*exp(k[1])))
    guess = log(length)*k[1]+k[2]
    Int(cld(bits, cld(bits, guess)))
end

extras(k) = [(
    a = x[2];
    b = heuristic_k(x[1]..., exp.(k.*100));
    b < firstindex(a) ? first(a)*2 : (b > lastindex(a) ? last(a)*2 : a[b]))/x[3]-1
    for x in data]
weights() = [(
    t = (x[1][1] == UInt ? 4 : 1) / sizeof(x[1][1]);
    sz = log10(x[1][3])+3;
    t*sz)
    for x in data]

function objective(k, full=false)
    ext = extras(k)
    ws = weights()

    miss_rate = sum(ws[ext .!= 0])/sum(ws)
    worst = maximum(ext)
    stdev = sqrt(sum(ws.*ext.^2)/sum(ws))
    meen = sum(ws.*ext)/sum(ws)
    fl = worst, stdev, meen, miss_rate
    w = 1, 3, 3, 1
    full ? fl : sum(fl .* w)
end

opt = optimize(objective, zeros(2))
    display(round.(objective(opt.minimizer, true).*1e5)./1e5)
    #(1.30975, 0.19471, 0.08808, 0.58667)
    k = exp.(opt.minimizer.*100)
    worst = data[sortperm(extras(opt.minimizer), rev=true)]
    worst_weights = weights()[sortperm(extras(opt.minimizer), rev=true)]
    x = sort(extras(opt.minimizer))
    histogram(x.+1)


##Check

struct RadixSort2ModifiedHeuristicAlg <: Algorithm
    delta::Int8
end

function sort!(v::AbstractVector, lo::Integer, hi::Integer, a::RadixSort2ModifiedHeuristicAlg, o::Ordering)

    T = Serializable(o, typeof(v))
    T === nothing && return sort!(v, lo, hi, MergeSort, o) # This has to be MergeSort rather
    # than defalg because defalg is unstable: https://github.com/JuliaLang/julia/issues/42713

    lo < hi || return v

    us, mn, mx = serialize(something(T), v, lo, hi, o)

    compression, bits, chunk_size = heuristic(mn, mx, hi-lo+1)

    us = compress!(us, compression)

    us = radix_sort!(us, similar(us, lo:hi), lo, hi, unsigned(bits),
        Val(UInt8(max(Int8(1), Int8(chunk_size)+Int8(a.delta)))))

    deserialize!(v, us, lo, hi, o, compression)
end

for T in sort(SortMark.UInts, by=sizeof)
    df = SortMark.make_df(
        [RadixSort2ModifiedHeuristicAlg(0),
            RadixSort2ModifiedHeuristicAlg(1),
            RadixSort2ModifiedHeuristicAlg(-1)],
        Types=[T], orders=[Base.Order.Forward],
        lens=SortMark.lengths(2, 1_000_000, 3),
        sources=Dict(:simple=>SortMark.sources[:simple]),
        seconds=nothing, samples=100)
    compute!(df)

    stat!(df, 1, 2)
    dfs = df[(!).(ismissing.(df.confint)), :]
    display(dfs[:, [:len, :Type, :pvalue, :point_estimate, :confint]])
    plot(log10.(dfs.len), first.(dfs.confint), label="+1",
        title="$(T): 95% ci of runtime ratio: normal/modified heuristic",
        xlabel="log10(length)", color=:red)
    plot!(log10.(dfs.len), last.(dfs.confint), label="+1", color=:red)

    stat!(df, 1, 3)
    dfs = df[(!).(ismissing.(df.confint)), :]
    println(T)
    display(dfs[:, [:len, :Type, :pvalue, :point_estimate, :confint]])
    plot!(log10.(dfs.len), first.(dfs.confint), label="-1", color=:blue)
    display(plot!(log10.(dfs.len), last.(dfs.confint), label="-1", color=:blue, ylims=[.7,1.3]))
end
