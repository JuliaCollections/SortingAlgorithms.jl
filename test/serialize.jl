@testset "serialize: quick" begin
    @test serialize(UInt8, Int8[3,1,2], 1, 3, Forward) == (UInt8[0x83, 0x81, 0x82], 0x81, 0x83)
    @test deserialize!(Vector{Int}(undef, 3), first(serialize(UInt, [3,1,2], 1, 3, Forward)), 1, 3, Forward, nothing) == [3,1,2]
end

@testset "serialize: full" begin
    #=for order in [Forward, Reverse], vals_t in vals
        ser, mn, mx = serialize(order, vals_t)
        @test mn == minimum(ser)
        @test mx == maximum(ser)
        deser = deserialize!(eltype(vals_t), order, copy(ser))
        @test all(deser .=== vals_t)
        sort!(ser)
        @test issorted(deserialize!(eltype(vals_t), order, ser), order=order)

        for a in vals_t
            @test serializable(order, eltype(vals_t))
            @test a === deserialize(eltype(vals_t), order, serialize(order, a))
            for b in vals_t
                a isa AbstractFloat && b isa AbstractFloat && isnan(a) && isnan(b) && continue # TODO revisit this departure from the standard
                @test lt(order, a, b) === (serialize(order, a) < serialize(order, b))
            end
        end
    end=#
end
