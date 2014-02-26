# SortingAlgorithms

[![Build Status](https://travis-ci.org/StefanKarpinski/SortingAlgorithms.jl.png)](https://travis-ci.org/StefanKarpinski/SortingAlgorithms.jl)

The SortingAlgorithms package provides three sorting algorithms that can be used with Julia's [standard sorting API](http://docs.julialang.org/en/latest/stdlib/sort/):

- [HeapSort] – an unstable, general purpose, in-place, O(n log n) comparison sort that works by heapify an array and repeatedly taking the maximal element from the heap.
- [TimSort] – a stable, general purpose, hybrid, O(n log n) comparison sort that adapts to different common patterns of partially ordered input data.
- [RadixSort] – a stable, special case, in-place, O(n) non-comparison sort that works by sorting data with fixed size, one digit at a time.

[HeapSort]:  http://en.wikipedia.org/wiki/Heapsort
[TimSort]:   http://en.wikipedia.org/wiki/Timsort
[RadixSort]: http://en.wikipedia.org/wiki/Radix_sort

## Usage

	julia> using SortingAlgorithms

	julia> words = map(chomp,[readlines(open("/usr/share/dict/words"))...])
	235886-element Array{ASCIIString,1}:
	 "A"
	 "a"
	 "aa"
	 "aal"
	 "aalii"
	 ⋮
	 "zythem"
	 "Zythia"
	 "zythum"
	 "Zyzomys"
	 "Zyzzogeton"

	julia> sort!(words, alg=TimSort)
	235886-element Array{ASCIIString,1}:
	 "A"
	 "Aani"
	 "Aaron"
	 "Aaronic"
	 "Aaronical"
	 ⋮
	 "zymotize"
	 "zymotoxic"
	 "zymurgy"
	 "zythem"
	 "zythum"

	julia> sort!(words, alg=TimSort, by=length)
	235886-element Array{ASCIIString,1}:
	 "A"
	 "B"
	 "C"
	 "D"
	 "E"
	 ⋮
	 "formaldehydesulphoxylate"
	 "pathologicopsychological"
	 "scientificophilosophical"
	 "tetraiodophenolphthalein"
	 "thyroparathyroidectomize"

	julia> sort!(words, alg=HeapSort)
	235886-element Array{ASCIIString,1}:
	 "A"
	 "Aani"
	 "Aaron"
	 "Aaronic"
	 "Aaronical"
	 ⋮
	 "zymotize"
	 "zymotoxic"
	 "zymurgy"
	 "zythem"
	 "zythum"

	julia> sort!(words, alg=HeapSort, by=length)
	235886-element Array{ASCIIString,1}:
	 "L"
	 "p"
	 "U"
	 "I"
	 "q"
	 ⋮
	 "pathologicopsychological"
	 "formaldehydesulphoxylate"
	 "scientificophilosophical"
	 "tetraiodophenolphthalein"
	 "thyroparathyroidectomize"

	julia> sort!(words, alg=RadixSort)
	ERROR: Radix sort only sorts bits types (got ASCIIString)
	 in error at error.jl:21
	 in sort! at /Users/stefan/.julia/SortingAlgorithms/src/SortingAlgorithms.jl:54
	 in sort! at sort.jl:328
	 in sort! at sort.jl:329

	julia> floats = randn(1000)
	1000-element Array{Float64,1}:
	  1.729
	  0.907196
	  0.461481
	 -0.204763
	 -0.16022
	  ⋮
	  0.700683
	 -0.236204
	 -2.15634
	 -0.316188
	 -0.171478

	julia> sort!(floats, alg=RadixSort)
	1000-element Array{Float64,1}:
	 -2.86255
	 -2.72041
	 -2.58234
	 -2.57259
	 -2.53046
	  ⋮
	  3.08307
	  3.12902
	  3.15075
	  3.20058
	  3.23942
