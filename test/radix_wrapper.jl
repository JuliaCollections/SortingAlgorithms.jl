using Test
using Random: shuffle

@testset "wrapped: quick" begin
    @test sort([3,1,2]) == sort([3,1,2], alg=RadixSort2)
    @test all(sort([NaN, -NaN]) .=== sort([NaN, -NaN], alg=RadixSort2))
    @test isequal(sort!([NaN], order=Base.Order.Reverse, alg=RadixSort2), [NaN])
    @test sort!([reinterpret(Char, 0xc1820000), 'a'], alg=RadixSort2) == ['a', reinterpret(Char, 0xc1820000)]
    @test sort!(["world", "hello"], alg=RadixSort2) == ["hello", "world"]
    x = vcat(Float16[NaN, -NaN], fill(0, 20))::Vector{Float16}
    @test all(sort(x; alg=MergeSort) .=== sort(x; alg=RadixSort2))
    @test_broken all(sort(x; alg=MergeSort) .=== sort(x))
end

@testset "wrapped: full" begin
    for oder in [Forward, Reverse]
        for vs in vals
            eltype(vs) == Float16 && continue
            vs = shuffle(vs)
            truth = sort(vs, order=oder)
            result = sort(vs, alg=RadixSort2, order=oder)
            @test typeof(truth) == typeof(result)
            if eltype(vs) <: AbstractFloat
                @test isequal(truth, result)
            else
                @test truth == result
            end
            @test all(truth .=== result)
            if !all(truth .=== result)
                println(truth, result)
            end
        end
    end
end
