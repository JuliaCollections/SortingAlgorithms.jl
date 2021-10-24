using Test
using OffsetArrays
using Base.Order

#TODO use soemone else's function
UIntN(::Val{8}) = UInt8
UIntN(::Val{16}) = UInt16
UIntN(::Val{32}) = UInt32
UIntN(::Val{64}) = UInt64
UIntN(::Val{128}) = UInt128

#Construct value lists
floats = [T[-π, -1.0, -1/π, 1/π, 1.0, π, -0.0, 0.0, Inf, -Inf, NaN, -NaN,
            prevfloat(T(0)), nextfloat(T(0)), prevfloat(T(Inf)), nextfloat(T(-Inf))]
    for T in [Float16, Float32, Float64]]

ints = [T[17, -T(17), 0, -one(T), 1, typemax(T), typemin(T), typemax(T)-1, typemin(T)+1]
    for T in [Int8, Int16, Int32, Int64, Int128, UInt8, UInt16, UInt32, UInt64, UInt128]]

vals = vcat(floats, ints)

#Add random vals
map(vals) do x
    append!(x, rand(eltype(x), 4))
    append!(x, reinterpret.(eltype(x), rand(UIntN(Val{sizeof(eltype(x))*8}()), 4)))
end

push!(vals, [true, false])

push!(vals, rand(Char, 20))

push!(vals, [7,3,9,1])
