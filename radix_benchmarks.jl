using SortingAlgorithms, SortMark

df = SortMark.make_df([QuickSort, MergeSort, RadixSort2])
expected = (df.Type .∉ [[UInt64, UInt128]]) .| (df.order .!= [Base.Order.Reverse]) .| (df.source_key .!= [:small_positive])
expected .&= (df.Type .∉ [[Float16,Char]])
df = df[expected, :]
compute!(df)
stat!(df)
