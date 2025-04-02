module SC
    export sc_polar_decoder, Properties, encode

    mutable struct Properties
        n::Int
        m::Int # log2(n)
        k::Int 

        frozen_bits::Array{UInt8, 1}

        function Properties(n::Int, k::Int, frozen_bits::Array{UInt8, 1})
            m::Int = Int(log2(n))

            @assert n > 0
            @assert k > 0
            
            p = new(n, m, k, frozen_bits)

            return p
        end
    end
    

        function f(a::Float64, b::Float64)::Float64
            return sign(a) * sign(b) * min(abs(a), abs(b))
        end
        
        function f(a::Array{Float64, 1},b::Array{Float64, 1})::Array{Float64, 1}
            return sign.(a) .* sign.(b) .* min.(abs.(a), abs.(b))
        end
        
        function g(a::Float64, b::Float64, u::UInt8)::Float64
            return b + (1 - 2 * Float64(u)) * a
        end
        
        function g(a::Array{Float64, 1}, b::Array{Float64, 1}, u::Array{UInt8, 1})::Array{Float64, 1}
            return b .+ (1 .- 2 * Float64.(u)) .* a
        end

        function encode(u::Array{UInt8, 1})
            if length(u) == 1
                return u
            else
                u2 = u[2:2:end]
                u1u2 = (u[1:2:end] + u2) .% 0x02
                return vcat(encode(u1u2), encode(u2))
            end
        end

        # attention: returns x_hat, not u_hat (run through polar encoder again to obtain u_hat)
        function sc_polar_decoder(LLRs::Array{Float64, 1}, frozen_bits::Array{UInt8, 1})
            N = length(LLRs)
            if N == 1
                # if frozen return 0 else return depending on LLR
                if frozen_bits[1] == 1
                    return [0x00]
                else
                    if LLRs[1] > 0
                        return [0x00]
                    else
                        return [0x01]
                    end
                end
            else
                # recursively decode
                half = div(N, 2)
                L1 = LLRs[1:half]
                L2 = LLRs[half+1:end]
                L_left = f(L1, L2)
                u1_hat = sc_polar_decoder(L_left, frozen_bits[1:half])
                L_rigth = g(L1, L2, u1_hat)
                u2_hat = sc_polar_decoder(L_rigth, frozen_bits[half+1:end])
                u1u2_hat = (u1_hat + u2_hat) .% 0x02
                return vcat(u1u2_hat, u2_hat)
            end
        end
end