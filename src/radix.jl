function radix_sort!(ps::AbstractVector{U}, ts::AbstractVector{U}, lo::Integer, hi::Integer, bits::Unsigned,
    ::Val{CHUNK_SIZE}) where {U <: Unsigned, CHUNK_SIZE}

    bits == 0 && return ps

    @inbounds begin

        MASK = UInt(1) << CHUNK_SIZE - 0x1
        counts = Vector{unsigned(typeof(hi-lo))}(undef, MASK+2)

        #ps0 = ps
        for shift in 0:CHUNK_SIZE:bits-1

            counts .= zero(eltype(counts))

            for k in lo:hi
                x = ps[k]
                idx = (x >> shift)&MASK + 2
                counts[idx] += one(eltype(counts))
            end

            sum = lo-1
            for i in 1:(MASK+1)
                sum += counts[i]
                counts[i] = sum
            end
            #accumulate!(+, counts, counts)

            for k in lo:hi # Is this iteration slower than it could be?
                x = ps[k]
                i = (x >> shift)&MASK + 1
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

function old_radix_sort!(ps::AbstractVector{U}, ts::AbstractVector{U}, length::Unsigned, bits::Unsigned,
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

#TODO clear thourough comment description
#TODO alternate counts as well as ps/ts
#TODO migrate regression testing target from old_radix_sort to radix_sort
#TODO deal with the fact that old_radix_sort is faster than radix_sort. Possibly use both?
# possibly drop floating point support?
# *Default arguments & comipler optimization with constants?
# *Use slower & leave a note that there is a faster availible for later performance optimizations.
# its only abt 5%
