using SortingAlgorithms, SortMark

df = SortMark.make_df([RadixSort2, QuickSort, MergeSort])
expected = (df.Type .∉ [[UInt64, UInt128]]) .| (df.order .!= [Base.Order.Reverse]) .| (df.source_key .!= [:small_positive])
df = df[expected, :]
compute!(df)
stat!(df, 1, 2)

x1 = df[(df.source_key .== :simple) .& (df.len .> 100) .& (df.Type .== Int) .&
    ((df.order .== [Base.Order.Forward]) .| (df.len .> 400)), :]
x1.seconds .= .05
compute!(x1)
stat!(x1, 1, 2)
println("always an improvement: ", maximum(first.(x1.confint)) < 1)
#Yay!

#Much still to do...
