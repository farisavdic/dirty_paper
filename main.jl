using Random, Statistics, Plots, LaTeXStrings, Dates

function polar_transform(u::Vector{Int})::Vector{Int}
    N = length(u)
    if N == 1
        return u
    end
    half = div(N, 2)
    u1 = [(u[i] ⊻ u[i + half]) for i in 1:half]
    u2 = u[half+1:end]
    return vcat(polar_transform(u1), polar_transform(u2))
end

function bpsk(x::Int)::Float64
    return 1 - 2 * x
end

function awgn(signal::Vector{Float64}, sigma::Float64)
    return signal .+ randn(length(signal)) .* sigma
end

function compute_llrs(y::Vector{Float64}, sigma::Float64)
    return (2 / sigma^2) .* y
end

function f(a::Float64, b::Float64)::Float64
    return sign(a) * sign(b) * min(abs(a), abs(b))
end

function g(a::Float64, b::Float64, u::Int)::Float64
    return b + (1 - 2 * u) * a
end

function sc_decode(llr::Vector{Float64}, frozen::Vector{Int}, index::Int = 1)::Vector{Int}
    N = length(llr)
    if N == 1
        return [frozen[index] == 1 ? 0 : (llr[1] >= 0 ? 0 : 1)]
    end
    half = div(N,2)
    llr_left = [f(llr[i], llr[i+half]) for i in 1:half]
    u_left = sc_decode(llr_left, frozen, index)
    llr_right = [g(llr[i], llr[i + half], u_left[i]) for i in 1:half]
    u_right = sc_decode(llr_right, frozen, index + half)
    return vcat(u_left, u_right)
end

# Simulation
function simulate_polar(N::Int, K::Int, sigma::Float64, trials::Int, mute::Bool, logs::Bool, log_path::String, simple_frozen::Bool)::Float64
    timestamp = now()

    # Define frozen bit positions (simple heuristic: freeze first N-K positions)
    if simple_frozen
        frozen = [i <= N-K ? 1 : 0 for i in 1:N]
    else
        # insert method to use proper unreliable channels
        return -1;
    end
    info_positions = findall(x -> x == 0, frozen)
    total_errors = 0
    frame_errors = 0

    for _ in 1:trials
        # Step 1: random info bits
        u_info = rand(0:1, K)

        # Step 2: Insert frozen bits
        u = zeros(Int, N)
        for (i, idx) in enumerate(info_positions)
            u[idx] = u_info[i]
        end

        # Step 3: Encode using polar transform
        x = polar_transform(u)

        # Step 4: BPSK modulation
        tx = [bpsk(bit) for bit in x]

        # Step 5: Pass through AWGN
        rx = awgn(tx, sigma)

        # Step 6: Compute LLRs
        llrs = compute_llrs(rx, sigma)

        # Step 7: Decode
        u_hat_full = sc_decode(llrs, frozen)

        # Step 8: Extract info bits
        u_hat = u_hat_full[info_positions]

        # Step 9: count errors
        if u_hat != u_info
            frame_errors += 1
        end
        total_errors += sum(u_hat .!= u_info)
    end
    fer = frame_errors / trials
    ber = total_errors / (K * trials)
    if logs
        data = "$(Dates.format(timestamp, "yyyy-mm-dd_HH:MM:SS"))\nN=$N\nK=$K\nsigma=$sigma\nsamples=$trials\nbit_errors=$total_errors\nber=$ber\nframe_errors=$frame_errors\nfer=$fer"
        open("$(log_path)$(Dates.format(timestamp, "yyyy-mm-dd_HH:MM:SS"))_log.txt","w") do f
            write(f, data)
        end
    end
    if !mute
        println("FER after $trials trials: $(round(fer, digits=5))")
        println("BER after $trials trials: $(round(ber, digits=5))")
    end
    return fer
end

#=
snr_dB = -2:0.01:2
snr_lin = zeros(length(snr_dB))
sigmas = zeros(length(snr_dB))
frame_error_rates = zeros(length(snr_dB))
for i in 1:length(snr_dB)
    snr_lin[i] = 10^(snr_dB[i]/10)
    sigmas[i] = sqrt(1 / snr_lin[i])
    frame_error_rates[i] = simulate_polar(1024,512,sigmas[i], 1000)
end
plot(snr_dB, frame_error_rates)
plot!(yscale=:log10, minorgrid=true)
xlims!(-2,2)
ylims!(1e-5, 1e+0)
xlabel!("SNR [dB]")
ylabel!("Frame Error Rate")
savefig("plot_1024.png")
=#

function simulation_sweep(N::Int, K::Int, samples::Int, lo_snr::Float64, hi_snr::Float64, step::Float64, mute::Bool, logs::Bool, log_path::String, simple_frozen::Bool)
    snr_dB = lo_snr:step:hi_snr
    snr_lin = zeros(length(snr_dB))
    sigmas = zeros(length(snr_dB))
    frame_error_rates = zeros(length(snr_dB))
    for i in 1:length(sigmas)
        snr_lin[i] = 10^(snr_dB[i]/10)
        sigmas[i] = sqrt(1 / snr_lin[i])
        frame_error_rates[i] = simulate_polar(N, K, sigmas[i], floor(Int, sqrt(i) * samples), true, true, "$log_path($N,$K)-sweep/", simple_frozen)
    end
end

function simulation_simple_frozen()
    lo_snr = -2.0
    hi_snr = 2.0
    step = 0.01
    log_path = "./data/simple_frozen/"
    sweeps = 10
    samples = 1000

    # (8,4)
    for i in 1:sweeps
        simulation_sweep(8,4, samples, lo_snr, hi_snr, step, true, true, log_path, true)
    end
    # (64,32)
    for i in 1:sweeps
        simulation_sweep(64,32, samples, lo_snr, hi_snr, step, true, true, log_path, true)
    end
    # (256, 128)
    for i in 1:sweeps
        simulation_sweep(256,128, samples, lo_snr, hi_snr, step, true, true, log_path, true)
    end
    # (1024, 512)
    for i in 1:sweeps
        simulation_sweep(1024,512, samples, lo_snr, hi_snr, step, true, true, log_path, true)
    end
end

simulation_simple_frozen()

