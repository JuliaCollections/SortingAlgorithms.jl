export QuadSort

struct QuadSortAlg <: Algorithm end

const QuadSort = maybe_optimize(QuadSortAlg())

# We remove all the branchless stuff. It actually doesn't do any good (we cannot do the clang-branchless optimization without
# introducing a lot of references, which makes everything much slower; and if we use the gcc-branchless workaround, the
# benchmarked performance is slightly worse).
# While branchless avoids branch mispredictions, we still have dependency chains; and Julia might sometimes decide to remove
# the branch anyway, but mimicking it with ifelse makes it harder for the compiler to do any optimization.
macro head_branchless_merge(arrayd, ptd, arrayl, ptl, arrayr, ptr, cmp)
    esc(quote
        let pl=$arrayl[$ptl], pr=$arrayr[$ptr]
            $arrayd[$ptd] = if $cmp(pl, pr) ≤ 0
                $ptl += 1
                pl
            else
                $ptr += 1
                pr
            end
        end
        $ptd += 1
    end)
end

macro tail_branchless_merge(arrayd, tpd, arrayl, tpl, arrayr, tpr, cmp)
    esc(quote
        let pl=$arrayl[$tpl], pr=$arrayr[$tpr]
            $arrayd[$tpd] = if $cmp(pl, pr) > 0
                $tpl -= 1
                pl
            else
                $tpr -= 1
                pr
            end
        end
        $tpd -= 1
    end)
end

macro getcmp(array, indexp1, fn, indexp2, cmp)
    esc(:(let a1=$array[$indexp1], a2=$array[$indexp2]
        ifelse($fn($cmp(a1, a2), 0), a1, a2) # here we keep it branchless, giving a minimal advantage over ternary
    end))
end

@inline function parity_merge_two!(array, array_index::Int, swap, swap_index::Int, cmp)
    @inbounds begin
        ptl = array_index; ptr = array_index + 2; pts = swap_index
        @head_branchless_merge(swap, pts, array, ptl, array, ptr, cmp)
        swap[pts] = @getcmp(array, ptl, ≤, ptr, cmp)
        ptl = array_index + 1; ptr = array_index + 3; pts = swap_index + 3
        @tail_branchless_merge(swap, pts, array, ptl, array, ptr, cmp)
        swap[pts] = @getcmp(array, ptl, >, ptr, cmp)
        nothing
    end
end

@inline function parity_merge_four!(array, array_index::Int, swap, swap_index::Int, cmp)
    @inbounds begin
        ptl = array_index; ptr = array_index + 4; pts = swap_index
        @head_branchless_merge(swap, pts, array, ptl, array, ptr, cmp)
        @head_branchless_merge(swap, pts, array, ptl, array, ptr, cmp)
        @head_branchless_merge(swap, pts, array, ptl, array, ptr, cmp)
        swap[pts] = @getcmp(array, ptl, ≤, ptr, cmp)
        ptl = array_index + 3; ptr = array_index + 7; pts = swap_index + 7
        @tail_branchless_merge(swap, pts, array, ptl, array, ptr, cmp)
        @tail_branchless_merge(swap, pts, array, ptl, array, ptr, cmp)
        @tail_branchless_merge(swap, pts, array, ptl, array, ptr, cmp)
        swap[pts] = @getcmp(array, ptl, >, ptr, cmp)
        nothing
    end
end

@inline function swap_branchless!(array, array_index::Int, cmp,
    x::Bool=@inbounds(cmp(array[array_index], array[array_index+1]) > 0))
    y = !x
    @inbounds array[array_index], array[array_index+1] = array[array_index+x], array[array_index+y]
    nothing
end

# the next seven functions are used for sorting 0 to 31 elements
function tiny_sort!(array, array_index::Int, nmemb::UInt, cmp)
    @inbounds if nmemb == 4
        # This is really just @inline quad_swap_four!(array, array_index, cmp)
        # However, for some weird reason, the profiler reports lots of GC and dynamic dispatch going on in either tiny_sort! or
        # quad_swap_four! unless we copy the code verbatim. Very strange and cannot be verified in the assembly code; however,
        # by inlining manually, it actually runs faster.
        swap_branchless!(array, array_index, cmp)
        swap_branchless!(array, array_index +2, cmp)
        array_index += 1

        a0 = array[array_index]
        a1 = array[array_index+1]
        if cmp(a0, a1) > 0
            array[array_index] = a1
            array[array_index+1] = a0
            array_index -= 1

            swap_branchless!(array, array_index, cmp)
            swap_branchless!(array, array_index +2, cmp)
            swap_branchless!(array, array_index +1, cmp)
        end
    elseif nmemb == 3
        swap_branchless!(array, array_index, cmp)
        swap_branchless!(array, array_index +1, cmp)
        swap_branchless!(array, array_index, cmp)
    elseif nmemb == 2
        swap_branchless!(array, array_index, cmp)
    end
end

# this function requires a minimum offset of 2 to work properly
function twice_unguarded_insert!(array, array_index::Int, offset::UInt, nmemb::UInt, cmp)
    @inbounds for i in asInt(offset):asInt(nmemb)-1
        end_ = array_index + i
        pta = end_
        cmp(array[pta -= 1], array[end_]) ≤ 0 && continue
        key = array[end_]
        if cmp(array[array_index+1], key) > 0
            top = i -1
            while true
                array[end_] = array[pta]; pta -= 1; end_ -= 1
                iszero(top -= 1) && break
            end
            array[end_] = key
            end_ -= 1
        else
            while true
                array[end_] = array[pta]; pta -= 1; end_ -= 1
                array[end_] = array[pta]; pta -= 1; end_ -= 1
                cmp(array[pta], key) > 0 || break
            end
            array[end_] = array[end_+1]
            array[end_+1] = key
        end
        swap_branchless!(array, end_, cmp)
    end
end

function quad_swap_four!(array, array_index::Int, cmp)
    @inbounds begin
        swap_branchless!(array, array_index, cmp)
        swap_branchless!(array, array_index +2, cmp)
        array_index += 1

        a0 = array[array_index]
        a1 = array[array_index+1]
        if cmp(a0, a1) > 0
            array[array_index] = a1
            array[array_index+1] = a0
            array_index -= 1

            swap_branchless!(array, array_index, cmp)
            swap_branchless!(array, array_index +2, cmp)
            swap_branchless!(array, array_index +1, cmp)
        end
        nothing
    end
end

function parity_swap_eight!(array, array_index::Int, swap, swap_index::Int, cmp)
    @inbounds begin
        swap_branchless!(array, array_index, cmp)
        swap_branchless!(array, array_index +2, cmp)
        swap_branchless!(array, array_index +4, cmp)
        swap_branchless!(array, array_index +6, cmp)

        cmp(array[array_index+1], array[array_index+2]) ≤ 0 &&
            cmp(array[array_index+3], array[array_index+4]) ≤ 0 &&
            cmp(array[array_index+5], array[array_index+6]) ≤ 0 && return
        parity_merge_two!(array, array_index, swap, swap_index, cmp)
        parity_merge_two!(array, array_index +4, swap, swap_index +4, cmp)
        parity_merge_four!(swap, swap_index, array, array_index, cmp)
    end
end

# left must be equal or one smaller than right
function parity_merge!(dest, dest_index::Int, from, from_index::Int, left::UInt, right::UInt, cmp)
    @inbounds begin
        ptl = from_index
        ptr = from_index + asInt(left)
        ptd = dest_index
        tpl = ptr -1
        tpr = tpl + asInt(right)
        tpd = dest_index + asInt(left + right) -1
        left < right && @head_branchless_merge(dest, ptd, from, ptl, from, ptr, cmp)
        @head_branchless_merge(dest, ptd, from, ptl, from, ptr, cmp)
        while !iszero(left -= 1)
            @head_branchless_merge(dest, ptd, from, ptl, from, ptr, cmp)
            @tail_branchless_merge(dest, tpd, from, tpl, from, tpr, cmp)
        end
        dest[tpd] = @getcmp(from, tpl, >, tpr, cmp)
        nothing
    end
end

function parity_swap_sixteen!(array, array_index::Int, swap, swap_index::Int, cmp)
    @inbounds begin
        quad_swap_four!(array, array_index, cmp)
        quad_swap_four!(array, array_index +4, cmp)
        quad_swap_four!(array, array_index +8, cmp)
        quad_swap_four!(array, array_index +12, cmp)
        cmp(array[array_index+3], array[array_index+4]) ≤ 0 &&
            cmp(array[array_index+7], array[array_index+8]) ≤ 0 &&
            cmp(array[array_index+11], array[array_index+12]) ≤ 0 && return
        parity_merge_four!(array, array_index, swap, swap_index, cmp)
        parity_merge_four!(array, array_index +8, swap, swap_index +8, cmp)
        parity_merge!(array, array_index, swap, swap_index, UInt(8), UInt(8), cmp)
    end
end

function tail_swap!(array, array_index::Int, swap, swap_index::Int, nmemb::UInt, cmp)
    @inbounds begin
        if nmemb < 5
            tiny_sort!(array, array_index, nmemb, cmp)
            return
        end
        if nmemb < 8
            quad_swap_four!(array, array_index, cmp)
            twice_unguarded_insert!(array, array_index, UInt(4), nmemb, cmp)
            return
        end
        if nmemb < 12
            parity_swap_eight!(array, array_index, swap, swap_index, cmp)
            twice_unguarded_insert!(array, array_index, UInt(8), nmemb, cmp)
            return
        end
        if 16 ≤ nmemb < 24
            parity_swap_sixteen!(array, array_index, swap, swap_index, cmp)
            twice_unguarded_insert!(array, array_index, UInt(16), nmemb, cmp)
            return
        end
        half1 = asInt(nmemb ÷ 2)
        quad1 = half1 ÷ 2
        quad2 = half1 - quad1

        half2 = asInt(nmemb - half1)
        quad3 = half2 ÷ 2
        quad4 = half2 - quad3

        pta = array_index
        tail_swap!(array, pta, swap, swap_index, asUInt(quad1), cmp); pta += quad1
        tail_swap!(array, pta, swap, swap_index, asUInt(quad2), cmp); pta += quad2
        tail_swap!(array, pta, swap, swap_index, asUInt(quad3), cmp); pta += quad3
        tail_swap!(array, pta, swap, swap_index, asUInt(quad4), cmp)

        cmp(array[array_index+quad1-1], array[array_index+quad1]) ≤ 0 &&
            cmp(array[array_index+half1-1], array[array_index+half1]) ≤ 0 &&
            cmp(array[pta-1], array[pta]) ≤ 0 && return

        parity_merge!(swap, swap_index, array, array_index, asUInt(quad1), asUInt(quad2), cmp)
        parity_merge!(swap, swap_index + half1, array, array_index + half1, asUInt(quad3), asUInt(quad4), cmp)
        parity_merge!(array, array_index, swap, swap_index, asUInt(half1), asUInt(half2), cmp)
    end
end

# the next three functions create sorted blocks of 32 elements
function quad_reversal!(array, pta::Int, ptz::Int)
    @inbounds begin
        loop = (ptz - pta) ÷ 2

        ptb = pta + loop
        pty = ptz - loop
        if iseven(loop)
            array[ptb], array[pty] = array[pty], array[ptb]
            ptb -= 1; pty += 1
            loop -= 1
        end
        loop ÷= 2
        while true
            array[pta], array[ptz] = array[ptz], array[pta]
            pta += 1; ptz -= 1
            array[ptb], array[pty] = array[pty], array[ptb]
            ptb -= 1; pty += 1

            iszero(loop) && break
            loop -= 1
        end
        nothing
    end
end

function quad_swap_merge!(array, array_index::Int, swap, swap_index::Int, cmp)
    parity_merge_two!(array, array_index, swap, swap_index, cmp)
    parity_merge_two!(array, array_index +4, swap, swap_index +4, cmp)
    parity_merge_four!(swap, swap_index, array, array_index, cmp)
end

function quad_swap!(array, array_index::Int, nmemb::UInt, cmp)
    @inbounds(@with_stackvec swap 32 eltype(array) begin
        pta = array_index
        count = nmemb ÷ 8
        while !iszero(count)
            count -= 1
            v1 = cmp(array[pta], array[pta+1]) > 0
            v2 = cmp(array[pta+2], array[pta+3]) > 0
            v3 = cmp(array[pta+4], array[pta+5]) > 0
            v4 = cmp(array[pta+6], array[pta+7]) > 0

            tmp = v1 + 2v2 + 4v3 + 8v4
            if tmp == 0
                cmp(array[pta+1], array[pta+2]) ≤ 0 &&
                    cmp(array[pta+3], array[pta+4]) ≤ 0 &&
                    cmp(array[pta+5], array[pta+6]) ≤ 0 &&
                    @goto ordered
                quad_swap_merge!(array, pta, swap, firstindex(swap), cmp)
            else
                if tmp == 15 && cmp(array[pta+1], array[pta+2]) > 0 &&
                    cmp(array[pta+3], array[pta+4]) > 0 &&
                    cmp(array[pta+5], array[pta+6]) > 0
                    pts = pta
                    @goto reversed
                end

                @label not_ordered
                swap_branchless!(array, pta, cmp, v1)
                swap_branchless!(array, pta +2, cmp, v2)
                swap_branchless!(array, pta +4, cmp, v3)
                swap_branchless!(array, pta +6, cmp, v4)

                quad_swap_merge!(array, pta, swap, firstindex(swap), cmp)
            end

            pta += 8
            continue

            @label ordered
            pta += 8
            if !iszero(count)
                count -= 1
                if (v1 = cmp(array[pta], array[pta+1]) > 0) | (v2 = cmp(array[pta+2], array[pta+3]) > 0) |
                    (v3 = cmp(array[pta+4], array[pta+5]) > 0) | (v4 = cmp(array[pta+6], array[pta+7]) > 0)
                    if v1 + v2 + v3 + v4 == 4 && cmp(array[pta+1], array[pta+2]) > 0 &&
                        cmp(array[pta+3], array[pta+4]) > 0 && cmp(array[pta+5], array[pta+6]) > 0
                        pts = pta
                        @goto reversed
                    end
                    @goto not_ordered
                end
                cmp(array[pta+1], array[pta+2]) ≤ 0 && cmp(array[pta+3], array[pta+4]) ≤ 0 &&
                    cmp(array[pta+5], array[pta+6]) ≤ 0 && @goto ordered
                quad_swap_merge!(array, pta, swap, firstindex(swap), cmp)
                pta += 8
                continue
            else
                count -= 1
            end
            break

            @label reversed
            pta += 8
            if !iszero(count)
                count -= 1
                if !((v1 = cmp(array[pta], array[pta+1]) ≤ 0) | (v2 = cmp(array[pta+2], array[pta+3]) ≤ 0) |
                    (v3 = cmp(array[pta+4], array[pta+5]) ≤ 0) | (v4 = cmp(array[pta+6], array[pta+7]) ≤ 0))
                    cmp(array[pta-1], array[pta]) > 0 && cmp(array[pta+1], array[pta+2]) > 0 &&
                        cmp(array[pta+3], array[pta+4]) > 0 && cmp(array[pta+5], array[pta+6]) > 0 &&
                        @goto reversed
                end
                quad_reversal!(array, pts, pta -1)
                v1 + v2 + v3 + v4 == 4 && cmp(array[pta+1], array[pta+2]) ≤ 0 &&
                    cmp(array[pta+3], array[pta+4]) ≤ 0 && cmp(array[pta+5], array[pta+6]) ≤ 0 &&
                    @goto ordered
                if v1 + v2 + v3 + v4 == 0 && cmp(array[pta+1], array[pta+2]) > 0 &&
                    cmp(array[pta+3], array[pta+4]) > 0 && cmp(array[pta+5], array[pta+6]) > 0
                    pts = pta
                    @goto reversed
                end

                swap_branchless!(array, pta, cmp, !v1)
                swap_branchless!(array, pta +2, cmp, !v2)
                swap_branchless!(array, pta +4, cmp, !v3)
                swap_branchless!(array, pta +6, cmp, !v4)

                if cmp(array[pta+1], array[pta+2]) > 0 ||
                    cmp(array[pta+3], array[pta+4]) > 0 ||
                    cmp(array[pta+5], array[pta+6]) > 0
                    quad_swap_merge!(array, pta, swap, firstindex(swap), cmp)
                end
                pta += 8
                continue
            else
                count -= 1
            end

            nmembrem = asInt(nmemb % 8) # subtracting -1 should give something negative
            doit = true
            for δ in nmembrem-1:-1:0
                if cmp(array[pta+δ-1], array[pta+δ]) ≤ 0
                    doit = false
                    break
                end
            end
            if doit
                quad_reversal!(array, pts, pta + nmembrem -1)
                pts == array_index && return true
                @goto reverse_end
            end
            quad_reversal!(array, pts, pta -1)
            break
        end
        tail_swap!(array, pta, swap, firstindex(swap), nmemb % 8, cmp)

        @label reverse_end

        pta = array_index
        count = nmemb ÷ 32
        while !iszero(count)
            count -= 1
            cmp(array[pta+7], array[pta+8]) ≤ 0 &&
                cmp(array[pta+15], array[pta+16]) ≤ 0 &&
                cmp(array[pta+23], array[pta+24]) ≤ 0 && continue
            parity_merge!(swap, firstindex(swap), array, pta, UInt(8), UInt(8), cmp)
            parity_merge!(swap, firstindex(swap) +16, array, pta +16, UInt(8), UInt(8), cmp)
            parity_merge!(array, pta, swap, firstindex(swap), UInt(16), UInt(16), cmp)

            pta += 32
        end
        nmemb % 32 > 8 && tail_merge!(array, pta, swap, firstindex(swap), UInt(32), nmemb % 32, UInt(8), cmp)
        return false
    end)
end

# quad merge support routines
function cross_merge!(dest, dest_index::Int, from, from_index::Int, left::UInt, right::UInt, cmp)
    @inbounds begin
        ptl = from_index
        ptr = from_index + asInt(left)
        tpl = ptr -1
        tpr = tpl + asInt(right)

        if left +1 ≥ right && right +1 ≥ left && left ≥ 32
            if cmp(from[ptl+15], from[ptr]) > 0 &&
                cmp(from[ptl], from[ptr+15]) ≤ 0 &&
                cmp(from[tpl], from[tpr-15]) > 0 &&
                cmp(from[tpl-15], from[tpr]) ≤ 0
                parity_merge!(dest, dest_index, from, from_index, left, right, cmp)
                return
            end
        end
        ptd = dest_index
        tpd = dest_index + asInt(left + right) -1

        while tpl - ptl > 8 && tpr - ptr > 8
            @label ptl8_ptr
            if cmp(from[ptl+7], from[ptr]) ≤ 0
                copyto!(dest, ptd, from, ptl, 8); ptd += 8; ptl += 8
                tpl - ptl > 8 && @goto ptl8_ptr
                break
            end
            @label ptl_ptr8
            if cmp(from[ptl], from[ptr+7]) > 0
                copyto!(dest, ptd, from, ptr, 8); ptd += 8; ptr += 8
                tpr - ptr > 8 && @goto ptl_ptr8
                break
            end
            @label tpl_tpr8
            if cmp(from[tpl], from[tpr-7]) ≤ 0
                tpd -= 8; tpr -= 8; copyto!(dest, tpd +1, from, tpr +1, 8)
                tpr - ptr > 8 && @goto tpl_tpr8
                break
            end
            @label tpl8_tpr
            if cmp(from[tpl-7], from[tpr]) > 0
                tpd -= 8; tpl -= 8; copyto!(dest, tpd +1, from, tpl +1, 8)
                tpl - ptl > 8 && @goto tpl8_tpr
            end
            loop = 8
            while true
                @head_branchless_merge(dest, ptd, from, ptl, from, ptr, cmp)
                @tail_branchless_merge(dest, tpd, from, tpl, from, tpr, cmp)
                iszero(loop -= 1) && break
            end
        end

        if cmp(from[tpl], from[tpr]) ≤ 0
            while ptl ≤ tpl
                @head_branchless_merge(dest, ptd, from, ptl, from, ptr, cmp)
            end
            while ptr ≤ tpr
                dest[ptd] = from[ptr]; ptr += 1; ptd += 1
            end
        else
            while ptr ≤ tpr
                @head_branchless_merge(dest, ptd, from, ptl, from, ptr, cmp)
            end
            while ptl ≤ tpl
                dest[ptd] = from[ptl]; ptl += 1; ptd += 1
            end
        end
        nothing
    end
end

# main memory: [A][B][C][D]
# swap memory: [A  B]       step 1
# swap memory: [A  B][C  D] step 2
# main memory: [A  B  C  D] step 3
function quad_merge_block!(array, array_index::Int, swap, swap_index::Int, block::UInt, cmp)
    @inbounds begin
        block_x_2 = 2block

        pt1 = array_index + asInt(block)
        pt2 = pt1 + asInt(block)
        pt3 = pt2 + asInt(block)

        tmp = (cmp(array[pt1-1], array[pt1]) ≤ 0) | 2(cmp(array[pt3-1], array[pt3]) ≤ 0)
        if tmp == 0
            cross_merge!(swap, swap_index, array, array_index, block, block, cmp)
            cross_merge!(swap, swap_index + asInt(block_x_2), array, pt2, block, block, cmp)
        elseif tmp == 1
            copyto!(swap, swap_index, array, array_index, block_x_2)
            cross_merge!(swap, swap_index + asInt(block_x_2), array, pt2, block, block, cmp)
        elseif tmp == 2
            cross_merge!(swap, swap_index, array, array_index, block, block, cmp)
            copyto!(swap, swap_index + asInt(block_x_2), array, pt2, block_x_2)
        elseif tmp == 3
            cmp(array[pt2-1], array[pt2]) ≤ 0 && return
            copyto!(swap, swap_index, array, array_index, 2block_x_2)
        end
        cross_merge!(array, array_index, swap, swap_index, block_x_2, block_x_2, cmp)
    end
end

function quad_merge!(array, array_index::Int, swap, swap_index::Int, swap_size::UInt, nmemb::UInt, block::UInt, cmp)
    pte = array_index + asInt(nmemb)
    block *= 4
    while block ≤ nmemb && block ≤ swap_size
        pta = array_index
        while true
            quad_merge_block!(array, pta, swap, swap_index, block ÷ 4, cmp)
            pta += asInt(block)
            pta + asInt(block) ≤ pte || break
        end
        tail_merge!(array, pta, swap, swap_index, swap_size, asUInt(pte - pta), block ÷ 4, cmp)
        block *= 4
    end
    tail_merge!(array, array_index, swap, swap_index, swap_size, nmemb, block ÷ 4, cmp)
    return block ÷ 2
end

function partial_forward_merge!(array, array_index::Int, swap, swap_index::Int, nmemb::UInt, block::UInt, cmp)
    @inbounds begin
        nmemb == block && return

        ptr = array_index + asInt(block)
        tpr = array_index + asInt(nmemb) -1
        cmp(array[ptr-1], array[ptr]) ≤ 0 && return

        copyto!(swap, swap_index, array, array_index, block)
        ptl = swap_index
        tpl = swap_index + asInt(block) -1
        while ptl < tpl -1 && ptr < tpr -1
            if cmp(swap[ptl], array[ptr+1]) > 0
                array[array_index] = array[ptr]
                array[array_index+1] = array[ptr+1]
                ptr += 2
                array_index += 2
            elseif cmp(swap[ptl+1], array[ptr]) ≤ 0
                array[array_index] = swap[ptl]
                array[array_index+1] = swap[ptl+1]
                ptl += 2
                array_index += 2
            else
                x = cmp(swap[ptl], array[ptr]) ≤ 0
                y = !x
                array[ptr+x] = array[ptr]
                ptr += 1
                array[array_index+y] = swap[ptl]
                ptl += 1
                array_index += 2
                @head_branchless_merge(array, array_index, swap, ptl, array, ptr, cmp)
            end
        end
        while ptl ≤ tpl && ptr ≤ tpr
            @head_branchless_merge(array, array_index, swap, ptl, array, ptr, cmp)
        end
        while ptl ≤ tpl
            array[array_index] = swap[ptl]
            ptl += 1
            array_index += 1
        end
        nothing
    end
end

function partial_backward_merge!(array, array_index::Int, swap, swap_index::Int, swap_size::UInt, nmemb::UInt, block::UInt,
    cmp)
    @inbounds begin
        nmemb == block && return

        tpl = array_index + asInt(block) -1
        tpa = array_index + asInt(nmemb) -1
        cmp(array[tpl], array[tpl+1]) ≤ 0 && return

        right = nmemb - block
        if nmemb ≤ swap_size && right ≥ 64
            cross_merge!(swap, swap_index, array, array_index, block, right, cmp)
            copyto!(array, array_index, swap, swap_index, nmemb)
            return
        end

        copyto!(swap, swap_index, array, array_index + asInt(block), right)
        tpr = swap_index + asInt(right) -1
        while tpl > array_index +16 && tpr > swap_index +16
            @label tpl_tpr16
            if cmp(array[tpl], swap[tpr-15]) ≤ 0
                loop = 16
                while true
                    array[tpa] = swap[tpr]
                    tpr -= 1
                    tpa -= 1
                    iszero(loop -= 1) && break
                end
                tpr > swap_index + 16 && @goto tpl_tpr16
                break
            end
            @label tpl16_tpr
            if cmp(array[tpl-15], swap[tpr]) > 0
                loop = 16
                while true
                    array[tpa] = array[tpl]
                    tpl -= 1
                    tpa -= 1
                    iszero(loop -= 1) && break
                end
                tpl > array_index +16 && @goto tpl16_tpr
                break
            end
            loop = 8
            while true
                if cmp(array[tpl], swap[tpr-1]) ≤ 0
                    array[tpa] = swap[tpr]
                    array[tpa-1] = swap[tpr-1]
                    tpr -= 2
                    tpa -= 2
                elseif cmp(array[tpl-1], swap[tpr]) > 0
                    array[tpa] = array[tpl]
                    array[tpa-1] = array[tpl-1]
                    tpl -= 2
                    tpa -= 2
                else
                    x = cmp(array[tpl], swap[tpr]) ≤ 0
                    y = !x
                    tpa -= 1
                    array[tpa+x] = swap[tpr]
                    tpr -= 1
                    array[tpa+y] = array[tpl]
                    tpl -= 1
                    tpa -= 1
                    @tail_branchless_merge(array, tpa, array, tpl, swap, tpr, cmp)
                end
                iszero(loop -= 1) && break
            end
        end
        while tpr > swap_index +1 && tpl > array_index +1
            if cmp(array[tpl], swap[tpr-1]) ≤ 0
                array[tpa] = swap[tpr]
                array[tpa-1] = swap[tpr-1]
                tpr -= 2
                tpa -= 2
            elseif cmp(array[tpl-1], swap[tpr]) > 0
                array[tpa] = array[tpl]
                array[tpa-1] = array[tpl-1]
                tpl -= 2
                tpa -= 2
            else
                x = cmp(array[tpl], swap[tpr]) ≤ 0
                y = !x
                tpa -= 1
                array[tpa+x] = swap[tpr]
                tpr -= 1
                array[tpa+y] = array[tpl]
                tpl -= 1
                tpa -= 1
                @tail_branchless_merge(array, tpa, array, tpl, swap, tpr, cmp)
            end
        end
        while tpr ≥ swap_index && tpl ≥ array_index
            @tail_branchless_merge(array, tpa, array, tpl, swap, tpr, cmp)
        end
        while tpr ≥ swap_index
            array[tpa] = swap[tpr]
            tpr -= 1
            tpa -= 1
        end
        nothing
    end
end

function tail_merge!(array, array_index::Int, swap, swap_index::Int, swap_size::UInt, nmemb::UInt, block::UInt, cmp)
    pte = array_index + asInt(nmemb)
    while block < nmemb && block ≤ swap_size
        pta = array_index
        while pta + asInt(block) < pte
            if pta + 2asInt(block) < pte
                partial_backward_merge!(array, pta, swap, swap_index, swap_size, 2block, block, cmp)
                pta += 2asInt(block)
                continue
            end
            partial_backward_merge!(array, pta, swap, swap_index, swap_size, asUInt(pte - pta), block, cmp)
            break
        end
        block *= 2
    end
    nothing
end

# the next four functions provide in-place rotate merge support
function trinity_rotation!(array, array_index::Int, swap, swap_index::Int, swap_size::UInt, nmemb::UInt, left::UInt)
    @inbounds begin
        right = nmemb - left
        if swap_size > 65536
            swap_size = UInt(65536)
        end
        if left < right
            if left ≤ swap_size
                copyto!(swap, swap_index, array, array_index, left)
                copyto!(array, array_index, array, array_index + asInt(left), right)
                copyto!(array, array_index + asInt(right), swap, swap_index, left)
            else
                pta = array_index
                ptb = pta + asInt(left)
                bridge = right - left
                if bridge ≤ swap_size && bridge > 3
                    ptc = pta + asInt(right)
                    ptd = ptc + asInt(left)
                    copyto!(swap, swap_index, array, ptb, bridge)
                    for _ in 1:left
                        array[ptc -= 1] = array[ptd -= 1]
                        array[ptd] = array[ptb -= 1]
                    end
                    copyto!(array, pta, swap, swap_index, bridge)
                else
                    ptc = ptb
                    ptd = ptc + asInt(right)
                    for _ in 1:left÷2
                        ptb -= 1
                        ptd -= 1
                        array[pta], array[ptb], array[ptc], array[ptd] = array[ptc], array[pta], array[ptd], array[ptb]
                        pta += 1
                        ptc += 1
                    end
                    for _ in 1:(ptd - ptc)÷2
                        ptd -= 1
                        array[pta], array[ptc], array[ptd] = array[ptc], array[ptd], array[pta]
                        ptc += 1
                        pta += 1
                    end
                    for _ in 1:(ptd - pta)÷2
                        ptd -= 1
                        array[pta], array[ptd] = array[ptd], array[pta]
                        pta += 1
                    end
                end
            end
        elseif right < left
            if right ≤ swap_size
                copyto!(swap, swap_index, array, array_index + asInt(left), right)
                copyto!(array, array_index + asInt(right), array, array_index, left)
                copyto!(array, array_index, swap, swap_index, right)
            else
                pta = array_index
                ptb = pta + asInt(left)
                bridge = left - right
                if bridge ≤ swap_size && bridge > 3
                    ptc = pta + asInt(right)
                    ptd = ptc + asInt(left)
                    copyto!(swap, swap_index, array, ptc, bridge)
                    for _ in 1:right
                        array[ptc] = array[pta]
                        array[pta] = array[ptb]
                        ptc += 1
                        pta += 1
                        ptb += 1
                    end
                    copyto!(array, ptd - asInt(bridge), swap, swap_index, bridge)
                else
                    ptc = ptb
                    ptd = ptc + asInt(right)
                    for _ in 1:right÷2
                        ptb -= 1
                        ptd -= 1
                        array[pta], array[ptb], array[ptc], array[ptd] = array[ptc], array[pta], array[ptd], array[ptb]
                        pta += 1
                        ptc += 1
                    end
                    for _ in 1:(ptb - pta)÷2
                        ptb -= 1
                        ptd -= 1
                        array[pta], array[ptb], array[ptd] = array[ptd], array[pta], array[ptb]
                        pta += 1
                    end
                    for _ in 1:(ptd - pta)÷2
                        ptd -= 1
                        array[pta], array[ptd] = array[ptd], array[pta]
                        pta += 1
                    end
                end
            end
        else
            pta = array_index
            ptb = pta + asInt(left)
            for _ in 1:left
                array[pta], array[ptb] = array[ptb], array[pta]
                pta += 1
                ptb += 1
            end
        end
        nothing
    end
end

function monobound_binary_first!(array, array_index::Int, value, top::UInt, cmp)
    @inbounds begin
        end_ = array_index + asInt(top)
        while top > 1
            mid = asInt(top ÷ 2)
            if cmp(value, array[end_-mid]) ≤ 0
                end_ -= mid
            end
            top -= mid
        end
        if cmp(value, array[end_-1]) ≤ 0
            end_ -= 1
        end
        return asUInt(end_ - array_index)
    end
end

function rotate_merge_block!(array, array_index::Int, swap, swap_index::Int, swap_size::UInt, lblock::UInt, right::UInt, cmp)
    @inbounds begin
        cmp(array[array_index+asInt(lblock)-1], array[array_index+asInt(lblock)]) ≤ 0 && return

        rblock = lblock ÷ 2
        rblocki = asInt(rblock)
        lblock -= rblock
        lblocki = asInt(lblock)

        left = monobound_binary_first!(array, array_index + lblocki + rblocki, array[array_index+lblocki], right, cmp)

        right -= left

        # [ lblock ] [ rblock ] [ left ] [ right ]

        if !iszero(left)
            if lblock + left ≤ swap_size
                copyto!(swap, swap_index, array, array_index, lblock)
                copyto!(swap, swap_index + lblocki, array, array_index + lblocki + rblocki, left)
                copyto!(array, array_index + lblocki + asInt(left), array, array_index + lblocki, rblock)

                cross_merge!(array, array_index, swap, swap_index, lblock, left, cmp)
            else
                trinity_rotation!(array, array_index + lblocki, swap, swap_index, swap_size, rblock + left, rblock)
                unbalanced = (2left < lblock) | (2lblock < left)
                if unbalanced && left ≤ swap_size
                    partial_backward_merge!(array, array_index, swap, swap_index, swap_size, lblock + left, lblock, cmp)
                elseif unbalanced && lblock ≤ swap_size
                    partial_forward_merge!(array, array_index, swap, swap_index, lblock + left, lblock, cmp)
                else
                    rotate_merge_block!(array, array_index, swap, swap_index, swap_size, lblock, left, cmp)
                end
            end
        end

        if !iszero(right)
            unbalanced = (2right < rblock) | (2rblock < right)
            if unbalanced && right ≤ swap_size || right + rblock ≤ swap_size
                partial_backward_merge!(array, array_index + lblocki + asInt(left), swap, swap_index, swap_size,
                    rblock + right, rblock, cmp)
            elseif unbalanced && rblock ≤ swap_size
                partial_forward_merge!(array, array_index + lblocki + asInt(left), swap, swap_index, rblock + right, rblock,
                    cmp)
            else
                rotate_merge_block!(array, array_index + lblocki + asInt(left), swap, swap_index, swap_size, rblock, right,
                    cmp)
            end
        end
        nothing
    end
end

function rotate_merge!(array, array_index::Int, swap, swap_index::Int, swap_size::UInt, nmemb::UInt, block::UInt, cmp)
    if nmemb ≤ 2block && nmemb - block ≤ swap_size # unsigned subtraction, ensures nmemb ≥ block
        partial_backward_merge!(array, array_index, swap, swap_index, swap_size, nmemb, block, cmp)
        return
    end

    pte = array_index + asInt(nmemb)
    while block < nmemb
        pta = array_index
        while pta + asInt(block) < pte
            if pta + 2asInt(block) < pte
                rotate_merge_block!(array, pta, swap, swap_index, swap_size, block, block, cmp)
                pta += 2asInt(block)
                continue
            end
            rotate_merge_block!(array, pta, swap, swap_index, swap_size, block, pte - pta - block, cmp)
            break
        end
        block *= 2
    end
end

function sort!(array::AbstractVector, lo::Int, hi::Int, ::QuadSortAlg, o::Ordering)
    len = UInt(hi - lo +1)
    cmp = let o=o; (x, y) -> lt(o, y, x) end

    if len < 32
        @with_stackvec(swap, 32, eltype(array),
            # we use a fixed size to make sure it is on the stack
            tail_swap!(array, lo, swap, firstindex(swap), len, cmp)
        )
    elseif !quad_swap!(array, lo, len, cmp)
        swap_size = len
        if swap_size ≤ 512
            @with_stackvec swapstack Int(swap_size) eltype(array) begin
                block = quad_merge!(array, lo, swapstack, firstindex(swapstack), swap_size, len, UInt(32), cmp)
                rotate_merge!(array, lo, swapstack, firstindex(swapstack), swap_size, len, block, cmp)
            end
        else
            local swap
            try
                swap = Vector{eltype(array)}(undef, swap_size)
            catch e
                if e isa OutOfMemoryError
                    @with_stackvec stack 512 eltype(array) begin
                        tail_merge!(array, lo, stack, firstindex(stack), UInt(32), len, UInt(32), cmp)
                        rotate_merge!(array, lo, stack, firstindex(stack), UInt(32), len, UInt(64), cmp)
                    end
                    return
                end
                rethrow()
            end
            block = quad_merge!(array, lo, swap, firstindex(swap), swap_size, len, UInt(32), cmp)
            rotate_merge!(array, lo, swap, firstindex(swap), swap_size, len, block, cmp)
        end
    end
    return array
end

function quadsort_swap!(array, array_index::Int, swap, swap_index::Int, swap_size::UInt, nmemb::UInt, cmp)
    if nmemb ≤ 96
        tail_swap!(array, array_index, swap, swap_index, nmemb, cmp)
    elseif !quad_swap!(array, array_index, nmemb, cmp)
        block = quad_merge!(array, array_index, swap, swap_index, swap_size, nmemb, UInt(32), cmp)
        rotate_merge!(array, array_index, swap, swap_index, swap_size, nmemb, block, cmp)
    end
end