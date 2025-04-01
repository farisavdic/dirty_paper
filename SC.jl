module SC
    export sc_polar_decoder, Properties, encode

    mutable struct Properties
        n::Int
        m::Int # log2(n)
        k::Int 

        frozen_bits::Dict{Int, UInt8}

        function Properties(n::Int, k::Int, frozen_bits::Dict{Int, UInt8})
            m::Int = Int(log2(n))

            @assert n > 0
            @assert k > 0
            
            p = new(n, m, k, frozen_bits)

            return p
        end

        function encode(u::Array{UInt8,1})
            if length(u) == 1
                x = u
            else
                u2 = u[2:2:end]
                u1u2 = ( u[1:2:end] + u2 ) .% 0x02
                x = [encode(u1u2); encode(u2)]
            end
        end

        function sc_polar_decoder(p::Properties, channel_LLRs::Array{Float64, 1})
            N = length(channel_LLRs)
            u_hat = zeros(UInt8, N)

            function recursive_decode(LLRs, start_idx, step)
                if step == 1
                    
                end
            end
        end
    end
end