using Printf, Dates
include("SCL.jl")

ps_5g(n) = 1 .+ filter(<(2^log2(n)), [0, 1, 2, 4, 8, 16, 32, 3, 5, 64, 9, 6, 17, 10, 18, 128, 12, 33, 65, 20, 256, 34, 24, 36, 7, 129, 66, 512, 11, 40, 68, 130, 19, 13, 48, 14, 72, 257, 21, 132, 35, 258, 26, 513, 80, 37, 25, 22, 136, 260, 264, 38, 514, 96, 67, 41, 144, 28, 69, 42, 516, 49, 74, 272, 160, 520, 288, 528, 192, 544, 70, 44, 131, 81, 50, 73, 15, 320, 133, 52, 23, 134, 384, 76, 137, 82, 56, 27, 97, 39, 259, 84, 138, 145, 261, 29, 43, 98, 515, 88, 140, 30, 146, 71, 262, 265, 161, 576, 45, 100, 640, 51, 148, 46, 75, 266, 273, 517, 104, 162, 53, 193, 152, 77, 164, 768, 268, 274, 518, 54, 83, 57, 521, 112, 135, 78, 289, 194, 85, 276, 522, 58, 168, 139, 99, 86, 60, 280, 89, 290, 529, 524, 196, 141, 101, 147, 176, 142, 530, 321, 31, 200, 90, 545, 292, 322, 532, 263, 149, 102, 105, 304, 296, 163, 92, 47, 267, 385, 546, 324, 208, 386, 150, 153, 165, 106, 55, 328, 536, 577, 548, 113, 154, 79, 269, 108, 578, 224, 166, 519, 552, 195, 270, 641, 523, 275, 580, 291, 59, 169, 560, 114, 277, 156, 87, 197, 116, 170, 61, 531, 525, 642, 281, 278, 526, 177, 293, 388, 91, 584, 769, 198, 172, 120, 201, 336, 62, 282, 143, 103, 178, 294, 93, 644, 202, 592, 323, 392, 297, 770, 107, 180, 151, 209, 284, 648, 94, 204, 298, 400, 608, 352, 325, 533, 155, 210, 305, 547, 300, 109, 184, 534, 537, 115, 167, 225, 326, 306, 772, 157, 656, 329, 110, 117, 212, 171, 776, 330, 226, 549, 538, 387, 308, 216, 416, 271, 279, 158, 337, 550, 672, 118, 332, 579, 540, 389, 173, 121, 553, 199, 784, 179, 228, 338, 312, 704, 390, 174, 554, 581, 393, 283, 122, 448, 353, 561, 203, 63, 340, 394, 527, 582, 556, 181, 295, 285, 232, 124, 205, 182, 643, 562, 286, 585, 299, 354, 211, 401, 185, 396, 344, 586, 645, 593, 535, 240, 206, 95, 327, 564, 800, 402, 356, 307, 301, 417, 213, 568, 832, 588, 186, 646, 404, 227, 896, 594, 418, 302, 649, 771, 360, 539, 111, 331, 214, 309, 188, 449, 217, 408, 609, 596, 551, 650, 229, 159, 420, 310, 541, 773, 610, 657, 333, 119, 600, 339, 218, 368, 652, 230, 391, 313, 450, 542, 334, 233, 555, 774, 175, 123, 658, 612, 341, 777, 220, 314, 424, 395, 673, 583, 355, 287, 183, 234, 125, 557, 660, 616, 342, 316, 241, 778, 563, 345, 452, 397, 403, 207, 674, 558, 785, 432, 357, 187, 236, 664, 624, 587, 780, 705, 126, 242, 565, 398, 346, 456, 358, 405, 303, 569, 244, 595, 189, 566, 676, 361, 706, 589, 215, 786, 647, 348, 419, 406, 464, 680, 801, 362, 590, 409, 570, 788, 597, 572, 219, 311, 708, 598, 601, 651, 421, 792, 802, 611, 602, 410, 231, 688, 653, 248, 369, 190, 364, 654, 659, 335, 480, 315, 221, 370, 613, 422, 425, 451, 614, 543, 235, 412, 343, 372, 775, 317, 222, 426, 453, 237, 559, 833, 804, 712, 834, 661, 808, 779, 617, 604, 433, 720, 816, 836, 347, 897, 243, 662, 454, 318, 675, 618, 898, 781, 376, 428, 665, 736, 567, 840, 625, 238, 359, 457, 399, 787, 591, 678, 434, 677, 349, 245, 458, 666, 620, 363, 127, 191, 782, 407, 436, 626, 571, 465, 681, 246, 707, 350, 599, 668, 790, 460, 249, 682, 573, 411, 803, 789, 709, 365, 440, 628, 689, 374, 423, 466, 793, 250, 371, 481, 574, 413, 603, 366, 468, 655, 900, 805, 615, 684, 710, 429, 794, 252, 373, 605, 848, 690, 713, 632, 482, 806, 427, 904, 414, 223, 663, 692, 835, 619, 472, 455, 796, 809, 714, 721, 837, 716, 864, 810, 606, 912, 722, 696, 377, 435, 817, 319, 621, 812, 484, 430, 838, 667, 488, 239, 378, 459, 622, 627, 437, 380, 818, 461, 496, 669, 679, 724, 841, 629, 351, 467, 438, 737, 251, 462, 442, 441, 469, 247, 683, 842, 738, 899, 670, 783, 849, 820, 728, 928, 791, 367, 901, 630, 685, 844, 633, 711, 253, 691, 824, 902, 686, 740, 850, 375, 444, 470, 483, 415, 485, 905, 795, 473, 634, 744, 852, 960, 865, 693, 797, 906, 715, 807, 474, 636, 694, 254, 717, 575, 913, 798, 811, 379, 697, 431, 607, 489, 866, 723, 486, 908, 718, 813, 476, 856, 839, 725, 698, 914, 752, 868, 819, 814, 439, 929, 490, 623, 671, 739, 916, 463, 843, 381, 497, 930, 821, 726, 961, 872, 492, 631, 729, 700, 443, 741, 845, 920, 382, 822, 851, 730, 498, 880, 742, 445, 471, 635, 932, 687, 903, 825, 500, 846, 745, 826, 732, 446, 962, 936, 475, 853, 867, 637, 907, 487, 695, 746, 828, 753, 854, 857, 504, 799, 255, 964, 909, 719, 477, 915, 638, 748, 944, 869, 491, 699, 754, 858, 478, 968, 383, 910, 815, 976, 870, 917, 727, 493, 873, 701, 931, 756, 860, 499, 731, 823, 922, 874, 918, 502, 933, 743, 760, 881, 494, 702, 921, 501, 876, 847, 992, 447, 733, 827, 934, 882, 937, 963, 747, 505, 855, 924, 734, 829, 965, 938, 884, 506, 749, 945, 966, 755, 859, 940, 830, 911, 871, 639, 888, 479, 946, 750, 969, 508, 861, 757, 970, 919, 875, 862, 758, 948, 977, 923, 972, 761, 877, 952, 495, 703, 935, 978, 883, 762, 503, 925, 878, 735, 993, 885, 939, 994, 980, 926, 764, 941, 967, 886, 831, 947, 507, 889, 984, 751, 942, 996, 971, 890, 509, 949, 973, 1000, 892, 950, 863, 759, 1008, 510, 979, 953, 763, 974, 954, 879, 981, 982, 927, 995, 765, 956, 887, 985, 997, 986, 943, 891, 998, 766, 511, 988, 1001, 951, 1002, 893, 975, 894, 1009, 955, 1004, 1010, 957, 983, 958, 987, 1012, 999, 1016, 767, 989, 1003, 990, 1005, 959, 1011, 1013, 895, 1006, 1014, 1017, 1018, 991, 1020, 1007, 1015, 1019, 1021, 1022, 1023])

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

function simulate_single_polar(p, sigma::Float64)
    bit_errors = 0
    frame_error = 0

    # generate k random data bits
    data = rand([0x00, 0x01], p.k)

    # combine data with frozen bits
    u = Array{UInt8}(undef, p.n)
    global j = 1
    for i=1:p.n
        if haskey(p.frozen_bits, i)
            u[i] = p.frozen_bits[i]
        else
            u[i] = data[j]
            global j += 1
        end
    end

    # polar transform
    #x = polar_transform(u)
    x = SCL.encode(u)

    # BPSK modulation and AWGN transmission
    y = (1 .- 2 * Float64.(x)) + randn(p.n) * sigma 

    # compute channel LLRs
    llrs = (2 / sigma ^ 2) .* y

    # polar decode
    #u_hat = sc_decode(llrs, frozen)
    @inbounds x_hat::Array{UInt8, 1} = SCL.scl_polar_decoder(p, llrs)

    # check result
    if x_hat != x
        bit_errors = sum(x .!= x_hat)
        frame_error = 1
        #@printf "c_hat != c, %i bit errors\n" bit_errors
    else
        #@printf "decoding successfull"
    end

    # return result vector
    return [frame_error, bit_errors]
end

function simulate_polar(p, snr_dB::Float64, samples::Int)
    bec = 0
    fec = 0
    ber = 0.0
    fer = 0.0

    sigma = sqrt(1 / (10^(snr_dB / 10)))
    for i in 1:samples
        tmp = simulate_single_polar(p, sigma)
        #println("SNR $snr_dB #$i")
        fec += tmp[1]
        bec += tmp[2]
    end
    ber = bec / (samples * k)
    fer = fec / samples
    res = [snr_dB, fer, ber, Float64(0.0), samples, fec, bec]
    @printf "SNR %.2f completed @ %s" snr_dB now()
    #println("SNR $snr_dB completed @ $(now())")
    return res
end

function read_test_case(filename)
    data = Dict()
    
    open(filename, "r") do file
        for line in eachline(file)
            key, value = split(line, ":", limit=2)
            key = strip(key)
            value = strip(value)

            if key == "n" || key == "k" || key == "nCRC"
                data[Symbol(key)] = parse(Int, value)
            elseif key == "frozen bits" || key == "crc polynomial"
                data[Symbol(replace(key, " " => "_"))] = parse.(Int, split(value, ","))
            elseif key == "design SNR"
                data[:design_SNR] = value == "None" ? nothing : parse(Float64, value)
            end
        end
    end
    
    return data
end

function run_simulations(test_directory_path::String, result_directory::String, L::Int, lo_snr::Float64, step::Float64)
    for file in readdir(test_directory_path)
        output_file = replace(file, "code"=>"results")
        out = "snr fer ber ber_c frames fec bec\n"
        config = read_test_case(string(test_directory_path, file))
        n = config[:n]
        k = config[:k]
        frz = config[:frozen_bits]
        CRCbits = config[:nCRC]
        CRC_hex = 0x01

        # pre-processing of frozen bit vector and CRC polynomial
        @assert n - k - CRCbits == sum(frz)
        num_frozen = Int(n - k - CRCbits)
        frozens = findall(vec(frz) .== 1)
        frozen_bits = Dict{Int, UInt8}()
        for i=1:num_frozen
            frozen_bits[frozens[i]] = 0x00
        end

        CRCpoly = Array{Bool}(undef, CRCbits+1)
        bitstr = bitstring(CRC_hex)[end-CRCbits+1:end]
        for j = 1:CRCbits
            if bitstr[j] == '1'
                CRCpoly[j] = true
            else
                CRCpoly[j] = false
            end
        end
        CRCpoly[CRCbits+1] = true

        p = SCL.Properties(n, k, L, frozen_bits, CRCbits, CRCpoly)

        target_FER = 1e-3
        min_FC = 100000
        min_FEC = 100
        num_inner_iter = 10000
        max_frames = Int(1e7 * min_FEC)

        snr = lo_snr
        while true
            frame_count = 0
            frame_error_count = 0
            frame_error_rate = 0.0
            bit_error_count = 0
            sigma = sqrt(1 / (10^(snr / 10)))
            while frame_count < min_FC || frame_error_count < min_FEC
                for _ in 1:num_inner_iter
                    tmp = simulate_single_polar(p, sigma)
                    frame_error_count += tmp[1]
                    bit_error_count += tmp[2]
                end
                frame_count += num_inner_iter
                frame_error_rate = frame_error_count / frame_count
            end
            if frame_error_rate <= target_FER
                break
            else
                bit_error_rate = bit_error_count / config[:k] * frame_count
                out *= @sprintf "%.6f %.6e %.6e %.6e %d %d %d\n" snr frame_error_rate bit_error_rate 0.0 frame_count frame_error_count bit_error_count
                @printf "SNR %.2f completed @ %s\n" snr Dates.format(now(), "HH:MM:SS")
                if frame_count >= max_frames
                    break
                else
                    snr += step
                end
            end
        end

        # test case finished
        println(out)
        #@printf out
        open(string(result_directory, output_file), "w") do f
            write(f, out)
        end
    end
end

run_simulations("prot/", "", 1, -1.0, 0.25)