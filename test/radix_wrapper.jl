using Test
using Random: shuffle

@testset "wrapped: quick" begin
    @test sort([3,1,2]) == sort([3,1,2], alg=RadixSort2)
    @test all(sort([NaN, -NaN]) .=== sort([NaN, -NaN], alg=RadixSort2))
    @test sort!([NaN], order=Base.Order.Reverse, alg=RadixSort2) == [NaN]
    #@test all(sort(Float16[NaN, -NaN]) .=== sort(Float16[NaN, -NaN], alg=RadixSort2))
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
