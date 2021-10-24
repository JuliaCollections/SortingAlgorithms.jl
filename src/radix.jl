function radix_sort!(ps::AbstractVector{U}, ts::AbstractVector{U}, length::Unsigned, bits::Unsigned,
    chunk_size::Unsigned) where {U <: Unsigned}

    bits == 0 && return ps

    @inbounds begin

        counts = Vector{typeof(length)}(undef, 2^chunk_size+1)
        mask = 2^chunk_size-1

        #ps0 = ps
        for shift in 0:chunk_size:bits-1

            counts .= zero(length)

            for x in ps
                idx = (x >> shift)&mask + 2
                counts[idx] += one(length)
            end

            sum = counts[1]
            for i in 2:2^chunk_size
                sum += counts[i]
                counts[i] = sum
            end
            #accumulate!(+, counts, counts)

            for x in ps
                i = (x >> shift)&mask + 1
                j = counts[i] += 1
                ts[j] = x
                #sortperm here
            end

            ps, ts = ts, ps

        end

        #if ps0 !== ps
        #    unsafe_copyto!(ps0, 1, ps, 1, length)
        #end

        ps#0
    end
end
