@testset "radix: quick" begin
    @test radix_sort!(UInt[3,1,2], Vector{UInt}(undef, 3), 1, 3, 0x2, Val(0x2)) == UInt[1,2,3]
end

using Random: shuffle
@testset "radix: full" begin
    for vs in vals
        if eltype(vs) <: Unsigned
            vs = shuffle(vs)
            truth = sort(vs)
            result = radix_sort!(vs, similar(vs), firstindex(vs), lastindex(vs), unsigned(sizeof(eltype(vs))*8), Val(0x9))
            @test typeof(truth) == typeof(result)
            @test truth == result
            @test all(truth .=== result)
        end
    end
end
