module SC
    export sc_polar_decoder, Properties, encode, bitrev

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

        function kron_power(A, n)
            if n == 1
                return A
            else
                return kron(A, kron_power(A, n-1))
            end
        end

        function bitrev(u)
            if length(u) == 2
                u_rev = u
            else
                u_rev = [bitrev(u[1:2:end]); bitrev(u[2:2:end])]
            end
            return u_rev
        end

        function encode(u::Array{UInt8, 1}, bit_reverse::Bool=false)
            if bit_reverse
                if length(u) == 1
                    return u
                else
                    u2 = u[2:2:end]
                    u1u2 = (u[1:2:end] + u2) .% 0x02
                    return vcat(encode(u1u2, true), encode(u2, true))
                end
            else
                G_2 = [1 0; 1 1]
                G = kron_power(G_2, log2(length(u)))
                return (transpose(transpose(u) * G) .% 2)
            end
        end

        # attention: returns x_hat, not u_hat (run through polar encoder again to obtain u_hat)
        function sc_polar_decoder(LLRs::Array{Float64, 1}, frozen_bits::Array{UInt8, 1}, bit_reverse::Bool=false)
            N = length(LLRs)

            if N == 1
                # not sure if this part works like this -> find way to get soft information!!!
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
                L_left = zeros(Float64, half)
                L_right = zeros(Float64, half)
                if bit_reverse
                    for i in 1:half
                        L_left[i] = f(LLRs[(2*i)-1], LLRs[2*i])
                    end
                    #L_left = SC.bitrev(L_left)
                    u1_hat = sc_polar_decoder(L_left, frozen_bits[1:half], true)
                    for i in 1:half
                        L_right[i] = g(LLRs[(2*i)-1], LLRs[2*i], u1_hat[i])
                    end
                    #L_right = SC.bitrev(L_right)
                    u2_hat = sc_polar_decoder(L_right, frozen_bits[half+1:end], true)
                    u1u2_hat = (u1_hat .+ u2_hat) .% 0x02
                    return vcat(u1u2_hat, u2_hat)
                else
                    for i in 1:half
                        L_left[i] = f(LLRs[i], LLRs[i+half])
                    end
                    u1_hat = sc_polar_decoder(L_left, frozen_bits[1:half])
                    for i in 1:half
                        L_right[i] = g(LLRs[i], LLRs[i+half], u1_hat[i])
                    end
                    u2_hat = sc_polar_decoder(L_right, frozen_bits[half+1:end])
                    u1u2_hat = (u1_hat .+ u2_hat) .% 0x02
                    return vcat(u1u2_hat, u2_hat)
                end
            end
        end

end