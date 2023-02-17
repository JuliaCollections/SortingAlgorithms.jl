using SortingAlgorithms
using Test
using StatsBase
using Random

a = rand(1:10000, 1000)

for alg in [TimSort, HeapSort, RadixSort, CombSort, BranchyPatternDefeatingQuicksort, BranchlessPatternDefeatingQuicksort]
    b = sort(a, alg=alg)
    @test issorted(b)
    ix = sortperm(a, alg=alg)
    b = a[ix]
    @test issorted(b)
    @test a[ix] == b

    b = sort(a, alg=alg, rev=true)
    @test issorted(b, rev=true)
    ix = sortperm(a, alg=alg, rev=true)
    b = a[ix]
    @test issorted(b, rev=true)
    @test a[ix] == b

    b = sort(a, alg=alg, by=x->1/x)
    @test issorted(b, by=x->1/x)
    ix = sortperm(a, alg=alg, by=x->1/x)
    b = a[ix]
    @test issorted(b, by=x->1/x)
    @test a[ix] == b

    c = copy(a)
    permute!(c, ix)
    @test c == b

    invpermute!(c, ix)
    @test c == a

    if alg != RadixSort  # RadixSort does not work with Lt orderings
        c = sort(a, alg=alg, lt=(>))
        @test b == c
    end

    c = sort(a, alg=alg, by=x->1/x)
    @test b == c
end

randnans(n) = reinterpret(Float64,[rand(UInt64)|0x7ff8000000000000 for i=1:n])

function randn_with_nans(n,p)
    v = randn(n)
    x = findall(rand(n).<p)
    v[x] = randnans(length(x))
    return v
end

Random.seed!(0xdeadbeef)

for n in [0:10..., 100, 101, 1000, 1001]
    r = 1:10
    v = rand(1:10,n)
    h = fit(Histogram, v, r)

    for ord in [Base.Order.Forward, Base.Order.Reverse]
        # insertion sort (stable) as reference
        pi = sortperm(v, alg=InsertionSort, order=ord)
        @test isperm(pi)
        si = v[pi]
        @test fit(Histogram, si, r) == h
        @test issorted(si, order=ord)
        @test all(issorted,[pi[si.==x] for x in r])
        c = copy(v)
        permute!(c, pi)
        @test c == si
        invpermute!(c, pi)
        @test c == v

        # stable algorithms
        for alg in [TimSort, RadixSort]
            p = sortperm(v, alg=alg, order=ord)
            @test p == pi
            s = copy(v)
            permute!(s, p)
            @test s == si
            invpermute!(s, p)
            @test s == v
        end

        # unstable algorithms
        for alg in [HeapSort, CombSort, BranchyPatternDefeatingQuicksort, BranchlessPatternDefeatingQuicksort]
            p = sortperm(v, alg=alg, order=ord)
            @test isperm(p)
            @test v[p] == si
            s = copy(v)
            permute!(s, p)
            @test s == si
            invpermute!(s, p)
            @test s == v
        end
    end

    v = randn_with_nans(n,0.1)
    for ord in [Base.Order.Forward, Base.Order.Reverse],
        alg in [TimSort, HeapSort, RadixSort, CombSort, BranchyPatternDefeatingQuicksort, BranchlessPatternDefeatingQuicksort]
        # test float sorting with NaNs
        s = sort(v, alg=alg, order=ord)
        @test issorted(s, order=ord)
        
        # This tests that NaNs (which compare equivalent) are treated stably 
        # even when the underlying algorithm is unstable. That it happens to
        # pass is not a part of the public API:
        @test reinterpret(UInt64, v[map(isnan, v)]) == reinterpret(UInt64, s[map(isnan, s)])

        # test float permutation with NaNs
        p = sortperm(v, alg=alg, order=ord)
        @test isperm(p)
        vp = v[p]
        @test isequal(vp,s)
        @test reinterpret(UInt64,vp) == reinterpret(UInt64,s)
    end
end

# additional tests to cover spacial cases of PdqSort
# test partial insertionsort, shuffle elements, partition_left
for v in [[1:1000;10], [1:500;500:-1:1], rand(Int,1000).%4]
    for alg in [BranchyPatternDefeatingQuicksort, BranchlessPatternDefeatingQuicksort]
        @test issorted(sort(v, alg=alg))
    end
end
# test fallback to HeapSort
let v = [1:500;500:-1:1]
    bad_allowed = 1
    SortingAlgorithms.pdqsort_loop!(v, 1, length(v), SortingAlgorithms.BranchyPatternDefeatingQuicksortAlg(), Base.Order.Forward, bad_allowed, nothing, nothing)
    @test issorted(v)
end
