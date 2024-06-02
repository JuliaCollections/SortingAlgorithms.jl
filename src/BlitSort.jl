export BlitSort

struct BlitSortAlg <: Algorithm end

const BlitSort = maybe_optimize(BlitSortAlg())

const BLIT_AUX = 512 # set to 0 for sqrt(n) cache size
const BLIT_OUT = 96 # should be smaller or equal to BLIT_AUX

function blit_analyze!(array, array_index::Int, swap, swap_index::Int, swap_size::UInt, nmemb::UInt, cmp::F) where {F}
    @inbounds begin
        half1 = nmemb ÷ 2
        quad1 = half1 ÷ 2
        quad2 = half1 - quad1
        half2 = nmemb - half1
        quad3 = half2 ÷ 2
        quad4 = half2 - quad3

        pta = array_index
        ptb = array_index + asInt(quad1)
        ptc = array_index + asInt(half1)
        ptd = array_index + asInt(half1 + quad3)

        astreaks = bstreaks = cstreaks = dstreaks = zero(UInt)
        abalance = bbalance = cbalance = dbalance = zero(UInt)

        cnt = nmemb
        while cnt > 132
            asum::UInt8 = bsum::UInt8 = csum::UInt8 = dsum::UInt8 = 0
            for _ in 32:-1:1
                asum += cmp(array[pta], array[pta+1]) > 0; pta += 1
                bsum += cmp(array[ptb], array[ptb+1]) > 0; ptb += 1
                csum += cmp(array[ptc], array[ptc+1]) > 0; ptc += 1
                dsum += cmp(array[ptd], array[ptd+1]) > 0; ptd += 1
            end
            abalance += asum; astreaks += asum = (asum == 0) | (asum == 32)
            bbalance += bsum; bstreaks += bsum = (bsum == 0) | (bsum == 32)
            cbalance += csum; cstreaks += csum = (csum == 0) | (csum == 32)
            dbalance += dsum; dstreaks += dsum = (dsum == 0) | (dsum == 32)

            if cnt > 516 && asum + bsum + csum + dsum == 0
                abalance += 48; pta += 96
                bbalance += 48; ptb += 96
                cbalance += 48; ptc += 96
                dbalance += 48; ptd += 96
                cnt -= 384
            end

            cnt -= 128
        end

        for _ in cnt:-4:UInt(8)
            abalance += cmp(array[pta], array[pta+1]) > 0; pta += 1
            bbalance += cmp(array[ptb], array[ptb+1]) > 0; ptb += 1
            cbalance += cmp(array[ptc], array[ptc+1]) > 0; ptc += 1
            dbalance += cmp(array[ptd], array[ptd+1]) > 0; ptd += 1
        end

        quad1 < quad2 && (bbalance += cmp(array[ptb], array[ptb+1]) > 0; ptb += 1)
        quad1 < quad3 && (cbalance += cmp(array[ptc], array[ptc+1]) > 0; ptc += 1)
        quad1 < quad4 && (dbalance += cmp(array[ptd], array[ptd+1]) > 0; ptd += 1)

        cnt = abalance + bbalance + cbalance + dbalance

        cnt == 0 && cmp(array[pta], array[pta+1]) ≤ 0 &&
            cmp(array[ptb], array[ptb+1]) ≤ 0 && cmp(array[ptc], array[ptc+1]) ≤ 0 && return

        abool = quad1 - abalance == 1
        bbool = quad2 - bbalance == 1
        cbool = quad3 - cbalance == 1
        dbool = quad4 - dbalance == 1

        if abool | bbool | cbool | dbool
            span1 = (abool && bbool) * (cmp(array[pta], array[pta+1]) > 0)
            span2 = (bbool && cbool) * (cmp(array[ptb], array[ptb+1]) > 0)
            span3 = (cbool && dbool) * (cmp(array[ptc], array[ptc+1]) > 0)

            tmp = span1 | span2 * 2 | span3 * 4
            if tmp == 1
                quad_reversal!(array, array_index, ptb)
                abalance = bbalance = 0
            elseif tmp == 2
                quad_reversal!(array, pta + 1, ptc)
                bbalance = cbalance = 0
            elseif tmp == 3
                quad_reversal!(array, array_index, ptc)
                abalance = bbalance = cbalance = 0
            elseif tmp == 4
                quad_reversal!(array, ptb + 1, ptd)
                cbalance = dbalance = 0
            elseif tmp == 5
                quad_reversal!(array, array_index, ptb)
                quad_reversal!(array, ptb + 1, ptd)
                abalance = bbalance = cbalance = dbalance = 0
            elseif tmp == 6
                quad_reversal!(array, pta + 1, ptd)
                bbalance = cbalance = dbalance = 0
            elseif tmp == 7
                quad_reversal!(array, array_index, ptd)
                return
            end

            abool && !iszero(abalance) && (quad_reversal!(array, array_index, pta); abalance = 0)
            bbool && !iszero(bbalance) && (quad_reversal!(array, pta + 1,     ptb); bbalance = 0)
            cbool && !iszero(cbalance) && (quad_reversal!(array, ptb + 1,     ptc); cbalance = 0)
            dbool && !iszero(dbalance) && (quad_reversal!(array, ptc + 1,     ptd); dbalance = 0)
        end

        cnt = nmemb ÷ 256 # more than 50% ordered
        abool = astreaks > cnt
        bbool = bstreaks > cnt
        cbool = cstreaks > cnt
        dbool = dstreaks > cnt

        tmp = abool + 2bbool + 4cbool + 8dbool
        if tmp == 0
            blit_partition!(array, array_index, swap, swap_index, swap_size, nmemb, cmp)
            return
        elseif tmp == 1
            iszero(abalance) || quadsort_swap!(array, array_index, swap, swap_index, swap_size, quad1, cmp)
            blit_partition!(array, pta + 1, swap, swap_index, swap_size, quad2 + half2, cmp)
        elseif tmp == 2
            blit_partition!(array, array_index, swap, swap_index, swap_size, quad1, cmp)
            iszero(bbalance) || quadsort_swap!(array, pta + 1, swap, swap_index, swap_size, quad2, cmp)
            blit_partition!(array, ptb + 1, swap, swap_index, swap_size, half2, cmp)
        elseif tmp == 3
            iszero(abalance) || quadsort_swap!(array, array_index, swap, swap_index, swap_size, quad1, cmp)
            iszero(bbalance) || quadsort_swap!(array, pta + 1, swap, swap_index, swap_size, quad2, cmp)
            blit_partition!(array, ptb + 1, swap, swap_index, swap_size, half2, cmp)
        elseif tmp == 4
            blit_partition!(array, array_index, swap, swap_index, swap_size, half1, cmp)
            iszero(cbalance) || quadsort_swap!(array, ptb + 1, swap, swap_index, swap_size, quad3, cmp)
            blit_partition!(array, ptc + 1, swap, swap_index, swap_size, quad4, cmp)
        elseif tmp == 8
            blit_partition!(array, array_index, swap, swap_index, swap_size, half1 + quad3, cmp)
            iszero(dbalance) || quadsort_swap!(array, ptc + 1, swap, swap_index, swap_size, quad4, cmp)
        elseif tmp == 9
            iszero(abalance) || quadsort_swap!(array, array_index, swap, swap_index, swap_size, quad1, cmp)
            blit_partition!(array, pta + 1, swap, swap_index, swap_size, quad2 + quad3, cmp)
            iszero(dbalance) || quadsort_swap!(array, ptc + 1, swap, swap_index, swap_size, quad4, cmp)
        elseif tmp == 12
            blit_partition!(array, array_index, swap, swap_index, swap_size, half1, cmp)
            iszero(cbalance) || quadsort_swap!(array, ptb + 1, swap, swap_index, swap_size, quad3, cmp)
            iszero(dbalance) || quadsort_swap!(array, ptc + 1, swap, swap_index, swap_size, quad4, cmp)
        else
            if abool
                iszero(abalance) || quadsort_swap!(array, array_index, swap, swap_index, swap_size, quad1, cmp)
            else
                blit_partition!(array, array_index, swap, swap_index, swap_size, quad1, cmp)
            end
            if bbool
                iszero(bbalance) || quadsort_swap!(array, pta + 1, swap, swap_index, swap_size, quad2, cmp)
            else
                blit_partition!(array, pta + 1, swap, swap_index, swap_size, quad2, cmp)
            end
            if cbool
                iszero(cbalance) || quadsort_swap!(array, ptb + 1, swap, swap_index, swap_size, quad3, cmp)
            else
                blit_partition!(array, ptb + 1, swap, swap_index, swap_size, quad3, cmp)
            end
            if dbool
                iszero(dbalance) || quadsort_swap!(array, ptc + 1, swap, swap_index, swap_size, quad4, cmp)
            else
                blit_partition!(array, ptc + 1, swap, swap_index, swap_size, quad4, cmp)
            end
        end

        if cmp(array[pta], array[pta+1]) ≤ 0
            if cmp(array[ptc], array[ptc+1]) ≤ 0
                cmp(array[ptb], array[ptb+1]) ≤ 0 && return
            else
                rotate_merge_block!(array, array_index + asInt(half1), swap, swap_index, swap_size, quad3, quad4, cmp)
            end
        else
            rotate_merge_block!(array, array_index, swap, swap_index, swap_size, quad1, quad2, cmp)
            cmp(array[ptc], array[ptc+1]) > 0 &&
                rotate_merge_block!(array, array_index + asInt(half1), swap, swap_index, swap_size, quad3, quad4, cmp)
        end
        rotate_merge_block!(array, array_index, swap, swap_index, swap_size, half1, half2, cmp)
    end
end

# The next 4 functions are used for pivot selection

function blit_binary_median(arraya, pta::Int, arrayb, ptb::Int, len::UInt, cmp::F) where {F}
    @inbounds begin
        while !iszero(len ÷= 2)
            leni = asInt(len)
            if cmp(arraya[pta+leni], arrayb[ptb+leni]) ≤ 0
                pta += leni
            else
                ptb += leni
            end
        end
        aa = arraya[pta]
        bb = arrayb[ptb]
        return ifelse(cmp(aa, bb) > 0, aa, bb)
    end
end

function blit_trim_four!(array, pta::Int, cmp::F) where {F}
    @inbounds begin
        x = cmp(array[pta], array[pta+1]) > 0
        array[pta], array[pta+1] = array[pta+x], array[pta+!x]
        pta += 2

        x = cmp(array[pta], array[pta+1]) > 0
        array[pta], array[pta+1] = array[pta+x], array[pta+!x]
        pta -= 2

        x = 2(cmp(array[pta], array[pta+2]) ≤ 0)
        array[pta+2] = array[pta+x]
        pta += 1
        x = 2(cmp(array[pta], array[pta+2]) > 0)
        array[pta] = array[pta+x]
    end
end

function blit_median_of_nine!(array, array_index::Int, swap, swap_index::Int, nmemb::UInt, cmp::F) where {F}
    @inbounds begin
        z = asInt(nmemb ÷ 9)

        pta = array_index

        for x in swap_index:swap_index+8
            swap[x] = array[pta]
            pta += z
        end

        blit_trim_four!(swap, swap_index, cmp)
        blit_trim_four!(swap, swap_index + 4, cmp)

        swap[swap_index]   = swap[swap_index+5]
        swap[swap_index+3] = swap[swap_index+8]

        blit_trim_four!(swap, swap_index, cmp)

        swap[swap_index] = swap[swap_index+6]

        x = cmp(swap[swap_index], swap[swap_index+1]) > 0
        y = cmp(swap[swap_index], swap[swap_index+2]) > 0
        z = cmp(swap[swap_index+1], swap[swap_index+2]) > 0

        return swap[swap_index + (x == y) + (y ⊻ z)]
    end
end

function blit_median_of_cbrt!(array, array_index::Int, swap, swap_index::Int, swap_size::UInt, nmemb::UInt, cmp::F) where {F}
    @inbounds begin
        cbrt = UInt(32) # TODO: figure out how to write this more efficiently using bsr
        while nmemb > cbrt * cbrt * cbrt && cbrt < swap_size
            cbrt *= 2
        end

        div = asInt(nmemb ÷ cbrt)

        pta = array_index
        pts = swap_index

        for cnt in 0:Core.bitcast(Int, cbrt)-1
            swap[pts+cnt] = array[pta]
            pta += div
        end
        pta = pts
        ptb = pts + asInt(cbrt ÷ 2)

        for cnt in cbrt÷8:-1:1
            blit_trim_four!(swap, pta, cmp)
            blit_trim_four!(swap, ptb, cmp)

            swap[pta] = swap[ptb+1]
            swap[pta+3] = swap[ptb+2]

            pta += 4
            ptb += 4
        end
        cbrt ÷= 4;

        quadsort_swap!(swap, pts, swap, pts + 2asInt(cbrt), cbrt, cbrt, cmp)
        quadsort_swap!(swap, pts + asInt(cbrt), swap, pts + 2asInt(cbrt), cbrt, cbrt, cmp)

        return cmp(swap[pts+2asInt(cbrt)-1], swap[pts]) ≤ 0,
            blit_binary_median(swap, pts, swap, pts + asInt(cbrt), cbrt, cmp)
    end
end

# As per suggestion by Marshall Lochbaum to improve generic data handling
function blit_reverse_partition!(array, array_index::Int, swap, swap_index::Int, swap_size::UInt, piv, nmemb::UInt, cmp::F) where {F}
    @inbounds begin
        if nmemb > swap_size
            h = nmemb ÷ 2;

            l = blit_reverse_partition!(array, array_index, swap, swap_index, swap_size, piv, h, cmp)
            r = blit_reverse_partition!(array, array_index + asInt(h), swap, swap_index, swap_size, piv, nmemb - h, cmp)

            trinity_rotation!(array, array_index + asInt(l), swap, swap_index, swap_size, h - l + r, h - l)

            return l + r
        end

        ptx = array_index
        pta = array_index
        pts = swap_index

        for _ in nmemb÷4:-1:1
            @unroll 4 begin
                if cmp(piv, array[ptx]) > 0
                    array[pta] = array[ptx]
                    pta += 1
                else
                    swap[pts] = array[ptx]
                    pts += 1
                end
                ptx += 1
            end
        end

        for _ in nmemb%4:-1:1
            if cmp(piv, array[ptx]) > 0
                array[pta] = array[ptx]
                pta += 1
            else
                swap[pts] = array[ptx]
                pts += 1
            end
            ptx += 1
        end
        m = pta - array_index

        _unsafe_copyto!(array, pta, swap, swap_index, nmemb - m)

        return asUInt(m)
    end
end

function blit_default_partition!(array, array_index::Int, swap, swap_index::Int, swap_size::UInt, piv, nmemb::UInt, cmp::F) where {F}
    @inbounds begin
        if nmemb > swap_size
            h = nmemb ÷ 2

            l = blit_default_partition!(array, array_index, swap, swap_index, swap_size, piv, h, cmp)
            r = blit_default_partition!(array, array_index + asInt(h), swap, swap_index, swap_size, piv, nmemb - h, cmp)

            trinity_rotation!(array, array_index + asInt(l), swap, swap_index, swap_size, h - l + r, h - l)

            return l + r
        end

        ptx = array_index
        pta = array_index
        pts = swap_index

        for _ in nmemb÷4:-1:one(UInt)
            @unroll 4 begin
                if cmp(array[ptx], piv) ≤ 0
                    array[pta] = array[ptx]
                    pta += 1
                else
                    swap[pts] = array[ptx]
                    pts += 1
                end
                ptx += 1
            end
        end

        for _ in nmemb%4:-1:one(UInt)
            if cmp(array[ptx], piv) ≤ 0
                array[pta] = array[ptx]
                pta += 1
            else
                swap[pts] = array[ptx]
                pts += 1
            end
            ptx += 1
        end
        m = pta - array_index

        _unsafe_copyto!(array, pta, swap, swap_index, nmemb - m)

        return asUInt(m)
    end
end

function blit_partition!(array, array_index::Int, swap, swap_index::Int, swap_size::UInt, nmemb::UInt, cmp::F) where {F}
    @inbounds begin
        a_size = zero(UInt)
        local max

        while true
            if nmemb ≤ 2048
                piv = blit_median_of_nine!(array, array_index, swap, swap_index, nmemb, cmp)
            else
                generic, piv = blit_median_of_cbrt!(array, array_index, swap, swap_index, swap_size, nmemb, cmp)

                if generic
                    quadsort_swap!(array, array_index, swap, swap_index, swap_size, nmemb, cmp)
                    return
                end
            end

            if !iszero(a_size) && cmp(max, piv) ≤ 0
                a_size = blit_reverse_partition!(array, array_index, swap, swap_index, swap_size, piv, nmemb, cmp)
                s_size = nmemb - a_size

                if s_size ≤ a_size÷16 || a_size <= BLIT_OUT
                    quadsort_swap!(array, array_index, swap, swap_index, swap_size, a_size, cmp)
                    return
                end
                nmemb = a_size
                a_size = zero(UInt)
                continue
            end

            a_size = blit_default_partition!(array, array_index, swap, swap_index, swap_size, piv, nmemb, cmp)
            s_size = nmemb - a_size

            if a_size ≤ s_size÷16 || s_size ≤ BLIT_OUT
                if iszero(s_size)
                    a_size = blit_reverse_partition!(array, array_index, swap, swap_index, swap_size, piv, a_size, cmp)
                    s_size = nmemb - a_size

                    if s_size ≤ a_size÷16 || a_size ≤ BLIT_OUT
                        return quadsort_swap!(array, array_index, swap, swap_index, swap_size, a_size, cmp)
                    end
                    nmemb = a_size
                    a_size = zero(UInt)
                    continue
                end
                quadsort_swap!(array, array_index + asInt(a_size), swap, swap_index, swap_size, s_size, cmp)
            else
                blit_partition!(array, array_index + asInt(a_size), swap, swap_index, swap_size, s_size, cmp)
            end

            if s_size ≤ a_size÷16 || a_size ≤ BLIT_OUT
                quadsort_swap!(array, array_index, swap, swap_index, swap_size, a_size, cmp)
                return
            end
            nmemb = a_size
            max = piv
        end
    end
end

function sort!(array::AbstractVector, lo::Int, hi::Int, ::BlitSortAlg, o::Ordering)
    len = UInt(hi - lo +1)
	if len ≤ 132
		sort!(array, lo, hi, QuadSortAlg(), o)
	else
        cmp = let o=o; (x, y) -> lt(o, y, x) end

        if !iszero(BLIT_AUX)
		    swap_size = UInt(BLIT_AUX)
        else
            swap_size = one(UInt) << 19

            while len ÷ swap_size < swap_size ÷ 128
                swap_size ÷= 4
            end
        end
        @with_stackvec swap Int(swap_size) eltype(array) begin
            blit_analyze!(array, lo, swap, firstindex(swap), swap_size, len, cmp)
        end
	end
    return array
end

blitsort_swap!(array, array_index::Int, swap, swap_index::Int, swap_size::UInt, nmemb::UInt, cmp::F) where {F} =
	(nmemb ≤ 132 ? quadsort_swap! : blit_analyze!)(array, array_index, swap, swap_index, swap_size, nmemb, cmp)