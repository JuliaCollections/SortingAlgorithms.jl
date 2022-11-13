# Sorting Algorithms

[![Build status](https://github.com/JuliaLang/SortingAlgorithms.jl/workflows/CI/badge.svg)](https://github.com/JuliaLang/SortingAlgorithms.jl/actions?query=workflow%3ACI+branch%3Amaster)
[![Coverage Status](https://coveralls.io/repos/JuliaLang/SortingAlgorithms.jl/badge.svg)](https://coveralls.io/r/JuliaLang/SortingAlgorithms.jl)
[![deps](https://juliahub.com/docs/SortingAlgorithms/deps.svg)](https://juliahub.com/ui/Packages/SortingAlgorithms/6dCmw?t=2)

The `SortingAlgorithms` package provides three sorting algorithms that can be used with Julia's [standard sorting API](https://docs.julialang.org/en/v1/base/sort/):

- [HeapSort] – an unstable, general purpose, in-place, O(n log n) comparison sort that works by heapifying an array and repeatedly taking the maximal element from the heap.
- [TimSort] – a stable, general purpose, hybrid, O(n log n) comparison sort that adapts to different common patterns of partially ordered input data.
- [CombSort] – an unstable, general purpose, in-place, O(n log n) comparison sort with O(n^2) pathological cases that is very efficient for primitieve types and orders due to SIMD instructions and instruction level paralellism on modern hardware.

[HeapSort]: https://en.wikipedia.org/wiki/Heapsort
[TimSort]:  https://en.wikipedia.org/wiki/Timsort
[CombSort]: https://en.wikipedia.org/wiki/Comb_sort

## Usage

```jl
	julia> using SortingAlgorithms

	julia> words = map(chomp,[readlines(open("/usr/share/dict/words"))...])
	235886-element Array{ASCIIString,1}:
	 "A"
	 "a"
	 "aa"
	 ⋮
	 "zythum"
	 "Zyzomys"
	 "Zyzzogeton"

	julia> sort!(words, alg=TimSort)
	235886-element Array{ASCIIString,1}:
	 "A"
	 "Aani"
	 "Aaron"
	 ⋮
	 "zymurgy"
	 "zythem"
	 "zythum"

	julia> sort!(words, alg=TimSort, by=length)
	235886-element Array{ASCIIString,1}:
	 "A"
	 "B"
	 "C"
	 ⋮
	 "scientificophilosophical"
	 "tetraiodophenolphthalein"
	 "thyroparathyroidectomize"

	julia> sort!(words, alg=HeapSort)
	235886-element Array{ASCIIString,1}:
	 "A"
	 "Aani"
	 "Aaron"
	 ⋮
	 "zymurgy"
	 "zythem"
	 "zythum"

	julia> sort!(words, alg=HeapSort, by=length)
	235886-element Array{ASCIIString,1}:
	 "L"
	 "p"
	 "U"
	 ⋮
	 "scientificophilosophical"
	 "tetraiodophenolphthalein"
	 "thyroparathyroidectomize"

	julia> sort!(randn(1000), alg=CombSort)
	1000-element Array{Float64,1}:
	 -2.86255
	 -2.72041
	 -2.58234
	  ⋮
	  3.15075
	  3.20058
	  3.23942
```

## Other packages that provide sorting algorithms

While SortingAlgorithms.jl is the most widely used sorting package in the Julia ecosystem, other packages are available:
- https://github.com/xiaodaigh/SortingLab.jl
- https://github.com/JeffreySarnoff/SortingNetworks.jl
- https://github.com/nlw0/ChipSort.jl
