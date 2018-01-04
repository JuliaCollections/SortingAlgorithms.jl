s = "abcdefghijklmno"
T = UInt128
skipbytes = 0
import SortingAlgorithms: StringRadixSort, StringRadixSortAlg
import Base.Perm

include("string_radix_sort.jl")
load_bits(UInt128,s)