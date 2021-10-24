using BenchmarkTools, Test, SortingAlgorithms

minmedtime(b) = (time(mean(b))+time(median(b)))/2
function evaluate(A, gen, V, L, bits, chunk, seconds=1, metric=minmedtime)

    o = @benchmark sort!(v; alg=$A) setup=(v=$gen()) evals=1 seconds=seconds

    n = @benchmark (eltype(v)).(radix_sort!($V.(v), Vector{$V}(undef, length(v)),
        $L(length(v)), UInt8($bits), $(UInt(unsigned(chunk))))) setup=(v=$gen()) evals=1 seconds=seconds

    metric(n) / metric(o)
end

trials = vcat([[
    (A, () -> rand(UInt64, 10_000), UInt64, UInt16, 64, 8),
    (A, () -> rand(1:1_000_000, 10_000), UInt32, UInt16, 20, 7),
    (A, () -> rand(1:100_000, 100), UInt32, UInt16, 17, 6),
    (A, () -> rand(UInt64, 1_000), UInt64, UInt16, 64, 4),
    (A, () -> rand(1:100, 100), UInt8, UInt8, 7, 4),
] for A in [QuickSort, RadixSort]]...)

standard = [0.71647, 0.37777, 0.70862, 1.1884, 0.51888, 1.10683, 0.91886, 0.04127, 0.47566, 0.03737]

result = [evaluate(t...) for t in trials]

pass = result .<= standard

all(pass) || println("rerunning $(count(.!pass)) regressions...")

result0 = copy(result)
result[.!pass] = [evaluate(t...) for t in trials[.!pass]]

pass = result .<= standard

let r = ceil.(min.(result, result0).*100_000)./100_000
    all(pass) || println("$(count(.!pass)) failed.
Lower the standards with: standard = $(max.(r, standard))")
end

@testset "regression" begin
    @test all(pass)
end
println("", "clearance: $(round.((1 .- (result ./ standard))*100*10)./10) %")
