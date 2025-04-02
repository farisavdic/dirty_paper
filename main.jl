include("SC.jl")
using Printf, Random
n = 8
k = 4
frozen_bits = [0x01, 0x01, 0x01, 0x00, 0x01, 0x00, 0x00, 0x00]
info_indices = [4, 6, 7, 8]
snr = 7.0
snr_lin = 10 ^ (snr / 10)
sigma = 1 / sqrt(snr_lin)

data = rand([0x00 0x01], k)
println("data: $data")

u = zeros(UInt8, n)
global c = 0
for i in info_indices
    u[i] = data[1 + c]
    global c += 1
end
println("u: $u")

x = SC.encode(u)
println("x: $x")

y = (1 .- 2 * Float64.(x)) + randn(n) * sigma
println("y: $y")

llrs = (2 / (sigma ^ 2)) * y
println("llrs: $llrs")
