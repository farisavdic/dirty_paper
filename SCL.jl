module SCL
	export scl_polar_decoder, Properties, encode
	
	## define properties of decoder in one struct
	mutable struct Properties
		n::Int	# code length
		m::Int	# log2(n)
		k::Int	# code dimension
		L::Int	# list size
		
		listhandover::Bool	# Are we part of a Multilevel Coding system and do we get an input list of LLR's?
		
		frozen_bits::Dict{Int,UInt8}
		CRCbits::Int
		CRCpoly::Array{Bool, 1}
		multipleCRC::Bool
		multipleCRC_interval::Int
		
		
		# Datastructures of SCL-Decoder
		inactivePathIndices::Array{Int, 1}
		activePath::Array{Bool, 1}
		arrays_L::Array{Array{Float64, 1}, 2}
		arrays_C::Array{Array{UInt8, 2}, 2}
		arrays_U::Array{UInt8, 2}
		pathIndexToArrayIndex::Array{Int, 2}
		inactiveArrayIndices::Array{Array{Int, 1}, 1}
		arrayReferenceCount::Array{Int, 2}
		pathMetric::Array{Float64, 2}
		tmp_PM::Array{Float64,2} # temporary PM storage so it doen't need to be reallocated all the time...
		unfrozen::Array{Int,1} # redundant info with frozen_bits, but simplifies computation...
		
		# only used when listhandover==true
		list_storage_idx::Array{Int,1} # index corresponding to data related to input list (for listhandover), update this index when cloning path
		pathMetric_handover::Array{Float64,1}
		
		function Properties(n::Int, k::Int, L::Int, frozen_bits::Dict{Int,UInt8}, CRCbits::Int=0, CRCpoly::Array{Bool,1}=[true], multipleCRC::Bool=false, multipleCRC_interval::Int=0, listhandover::Bool=false)
			m::Int = Int(log2(n))
			
			@assert L > 0
			@assert n > 0
			@assert k > 0
			@assert CRCbits >= 0
			
			if multipleCRC == false
				@assert k + CRCbits <= n
				@assert n == k + CRCbits + length(frozen_bits)
			else
				num_CRCs = ceil(k/multipleCRC_interval)
				@assert k + num_CRCs*CRCbits < n
				@assert n == k + length(frozen_bits) + num_CRCs*CRCbits
			end
			
			#init datastructures
			inactivePathIndices		= collect(1:L)
			activePath				= falses(L)
			arrays_L				= Array{Array{Float64, 1}, 2}(undef, m+1, L) # 2-D array of arrays
			arrays_C				= Array{Array{UInt8, 2}, 2}(undef, m+1, L)   # 2-D array of arrays
			arrays_U				= Array{UInt8, 2}(undef, n, L) # Additional Array for tracking the dataword U (including frozen bits)
			pathIndexToArrayIndex	= zeros(Int, m+1, L)
			inactiveArrayIndices	= Array{Array{Int, 1}, 1}(undef, m+1)
			arrayReferenceCount		= zeros(Int, m+1, L)
			pathMetric				= zeros(L, n)
			tmp_PM					= zeros(L, 2)
			
			@simd for lambda = 1:m+1
				@simd for s = 1:L
					arrays_L[lambda,s] = zeros(2^(m-lambda+1))	# float
					arrays_C[lambda,s] = Array{UInt8, 2}(undef, 2^(m-lambda+1), 2)
				end
				inactiveArrayIndices[lambda] = collect(1:L)
			end	
			
			# call constructor
			p = new(n, m, k, L, listhandover, frozen_bits, CRCbits, CRCpoly, multipleCRC, multipleCRC_interval, inactivePathIndices, activePath, arrays_L, arrays_C, arrays_U, pathIndexToArrayIndex, inactiveArrayIndices, arrayReferenceCount, pathMetric, tmp_PM)
			
			if listhandover # init index storage
				p.list_storage_idx = collect(1:L)
			end
			
			return p			
		end
	end

    # encode polar code without bit reversal
    function encode_nonBR(u)
        G2 = [1 0; 1 1]
        G = G2
        for i=1:log2(length(u))-1
            G = kron(G,G2)
        end
        
        x = (u*G) .% 2
        return x
    end

    # polar encoding including bit reversal
    function encode(u::Array{UInt8,1}) # following Alg. 1 in [Pf14] #TODO: optimize
        if length(u) == 1
            x = u
        else
            u2 = u[2:2:end]
            u1u2 = ( u[1:2:end] + u2 ) .% 0x02
            x = [encode(u1u2); encode(u2)]
        end
    end


    # bit reversal
    function bitrev(u)
        if length(u) == 2
            u_rev = u
        else
            u_rev = [bitrev(u[1:2:end]); bitrev(u[2:2:end])]
        end
        return u_rev
    end


	## calculates CRC of bits (Array of length n with contents 0x00 or 0x01 (UInt8))
	## CRC Polynomial as Bool, e.g. for size=8 and polynomial x^9 + x^8 + x^2 + x + 1: poly=[true true false false false false false true true true]
	function calc_CRC(bits::Array{UInt8}, size::UInt, poly::Array{Bool})
		b = copy(bits)
		append!(b, zeros(UInt8, size))
		
		for i = 1:length(bits)
			if popfirst!(b) == 0x01 #leading bit is one, xor the rest with polynomial
				@simd for j = 1:size
					b[j] = xor(poly[j+1], b[j])
				end
			end
		end
		
		return b
	end
	
	# [TaVa15, Alg. 6]
	function assignInitialPath(p::Properties)
		l = pop!(p.inactivePathIndices)
		p.activePath[l]	= true
		for lambda = 1:p.m+1
			s = pop!(p.inactiveArrayIndices[lambda])
			p.pathIndexToArrayIndex[lambda, l] = s
			p.arrayReferenceCount[lambda, s] = 1
		end
		
		return l
	end
	
	# [TaVa15, Alg. 9]
	function getArray_L(lambda::Int, l::Int, p::Properties)
		s = p.pathIndexToArrayIndex[lambda, l]
		
		if p.arrayReferenceCount[lambda, s] == 1
			s_prime = s
		else
			s_prime = pop!(p.inactiveArrayIndices[lambda])
			## copy 
			p.arrays_L[lambda, s_prime] = copy(p.arrays_L[lambda, s])
			p.arrays_C[lambda, s_prime] = copy(p.arrays_C[lambda, s])
			p.arrayReferenceCount[lambda, s] -= 1
			p.arrayReferenceCount[lambda, s_prime] += 1
			p.pathIndexToArrayIndex[lambda, l] = s_prime
		end
		
		return p.arrays_L[lambda, s_prime]
	end
	
	# [TaVa15, Alg. 9]
	function getArray_C(lambda::Int, l::Int, p::Properties)
		s = p.pathIndexToArrayIndex[lambda, l]
		
		if p.arrayReferenceCount[lambda, s] == 1
			s_prime = s
		else
			s_prime = pop!(p.inactiveArrayIndices[lambda])
			## copy 
			p.arrays_L[lambda, s_prime] = copy(p.arrays_L[lambda, s])
			p.arrays_C[lambda, s_prime] = copy(p.arrays_C[lambda, s])
			p.arrayReferenceCount[lambda, s] -= 1
			p.arrayReferenceCount[lambda, s_prime] += 1
			p.pathIndexToArrayIndex[lambda, l] = s_prime
		end
		
		return p.arrays_C[lambda, s_prime]
	end
	
	# [TaVa15, Alg. 10]
	# recursively calculates values of P arrays
	function recursively_calc_L(lambda::Int, phi::Int, m::Int, L::Int, p::Properties) # lambda from 1 to m+1
		if lambda == 1
			return # nothing more to do, end of recursion
		end
		
		Ψ::Int = phi >> 1 #Int(floor(phi / 2))
		if phi % 2 == 0
			recursively_calc_L(lambda-1, Ψ, m, L, p) # recurse (values haven't been computed yet)
		end
		
		# main calculation
		for l=1:L
			if p.activePath[l] # only continue if path is active
				Ln = getArray_L(lambda, l, p)
				L2 = getArray_L(lambda-1, l, p)
				C = getArray_C(lambda, l, p)
				
				for beta = 0:2^(m-lambda+1)-1
					if phi % 2 == 0 # case 1, equation (4) [TaVa15]
						Ln[beta+1] = f_minus(L2[2*beta+1], L2[2*beta+1+1])
					else # case 2, equation (5) [TaVa15]
						u = C[beta+1, 1]
						Ln[beta+1] = f_plus(L2[2*beta+1], L2[2*beta+1+1], u)
					end
				end
			end
		end
	end
	
	# recursively update C
	# [TaVa15, Alg. 11]
	function recursively_update_C(lambda::Int, phi::Int, m::Int, L::Int, p::Properties)
		if phi%2 != 1
			error("ERROR! recursively_update_C called with even Phi!")
			return
		end
		
		Ψ::Int = phi >> 1 # Int(floor(phi / 2))
		for l=1:L
			if p.activePath[l]
				C = getArray_C(lambda, l, p)
				C2 = getArray_C(lambda-1, l, p)
				@simd for beta = 0:2^(m-lambda+1) - 1
					C2[2*beta+1, Ψ%2+1] = ( C[beta+1, 1] + C[beta+1, 2] ) % 0x02
					C2[2*beta+1+1, Ψ%2+1] = C[beta+1, 2]
				end
			end
		end
		
		if Ψ % 2 == 1
			recursively_update_C(lambda-1, Ψ, m, L, p)
		end
	end
	
	# clone a path and return index of new copy
	# [TaVa15, Alg. 7]
	function clonePath(l::Int, m::Int, phi::Int, p::Properties)
		lprime = pop!(p.inactivePathIndices)
		p.activePath[lprime] = true
		
		# set references
		@simd for lambda = 1:m+1
			s = p.pathIndexToArrayIndex[lambda, l]
			p.pathIndexToArrayIndex[lambda, lprime] = s
			p.arrayReferenceCount[lambda,s] += 1
		end
		
		# set arrays_U
		p.arrays_U[1:phi+1, lprime] = p.arrays_U[1:phi+1, l]
		
		# update list storage indizes
		if p.listhandover
			@assert p.list_storage_idx[lprime] == 0 #check if not overwriting existing storage
			p.list_storage_idx[lprime] = p.list_storage_idx[l]
		end
		
		return lprime
	end
	
	
	# kill a Path
	# [TaVa15, Alg. 8]
	function killPath(l::Int, m::Int, p::Properties)
		p.activePath[l] = false
		push!(p.inactivePathIndices, l)
		
		for lambda = 1:m+1
			s = p.pathIndexToArrayIndex[lambda, l]
			p.arrayReferenceCount[lambda, s] -= 1
			if p.arrayReferenceCount[lambda, s] == 0
				push!(p.inactiveArrayIndices[lambda], s)
			end
		end
		
		# update list storage indizes
		if p.listhandover
			p.list_storage_idx[l] = 0 # set storage index to 0 (=not used)
		end
	end
	
	
	# [StPaBu15]
	function calc_PM(PMprime::Float64, L::Float64, u::UInt8)
		#PM = PMprime + log(1 + e^(-(1-2*u) * L) )
		if u == (1 - sign(L))/2 #approximation
			PM = PMprime
		else
			PM = PMprime + abs(L)
		end
	end
	
	# [StPaBu15]
	function f_minus(a::Float64, b::Float64)
		#f = log( (exp(a+b) + 1) / (exp(a) + exp(b)) )
		f = sign(a) * sign(b) * min( abs(a), abs(b) ) #approximation
	end
	
	function f_plus(a::Float64, b::Float64, u::UInt8)
		f = (-1)^u * a + b
	end
	
	
	# [TaVa15, Alg. 13], combined with [StPaBu15]
	function continuePaths_UnfrozenBit(phi::Int, m::Int, L::Int, p::Properties)
		
		if 2*sum(p.activePath) <= L #just copy all active paths
			acopy = copy(p.activePath)
			for l = 1:L
				if acopy[l]
					C = getArray_C(m+1, l, p)
					Ln = getArray_L(m+1, l, p)
					
					C[1, phi%2+1] = 0x00
					p.arrays_U[phi+1, l] = 0x00
					if phi > 0
						p.pathMetric[l, phi+1] = calc_PM(p.pathMetric[l, phi], Ln[1], 0x00)
					else	
						if p.listhandover
							p.pathMetric[l, phi+1] = calc_PM(p.pathMetric_handover[l], Ln[1], 0x00)
						else
							p.pathMetric[l, phi+1] = calc_PM(0.0, Ln[1], 0x00)
						end
					end
					
					lprime = clonePath(l, m, phi, p)
					C = getArray_C(m+1, lprime, p)
					C[1, phi%2+1] = 0x01
					if phi > 0
						p.pathMetric[lprime, phi+1] = calc_PM(p.pathMetric[l, phi], Ln[1], 0x01)
					else
						if p.listhandover
							p.pathMetric[l, phi+1] = calc_PM(p.pathMetric_handover[l], Ln[1], 0x01)
						else
							p.pathMetric[lprime, phi+1] = calc_PM(0.0, Ln[1], 0x01)
						end
					end
					p.arrays_U[phi+1, lprime] = 0x01
				end
			end
			
		else # some paths have to be removed first
			# calculate PathMetric for all possible paths
			@simd for l = 1:L
				if p.activePath[l]
					Ln = getArray_L(m+1, l, p)
					@simd for s = 0:1
						if phi > 0
							p.tmp_PM[l, s+1] = calc_PM(p.pathMetric[l, phi], Ln[1], UInt8(s))
						else
							if p.listhandover
								p.tmp_PM[l, s+1] = calc_PM(p.pathMetric_handover[l], Ln[1], UInt8(s))
							else
								p.tmp_PM[l, s+1] = calc_PM(0.0, Ln[1], UInt8(s))
							end
						end
					end
				end
			end
			
			vPM = vec(p.tmp_PM)
			keepidx = sortperm(vPM)[1:L]
			keep = falses(L,2);
			for l = 1:L
				if keepidx[l] <= L
					keep[keepidx[l],1] = true
				else
					keep[keepidx[l]-L,2] = true
				end
			end
			
			# kill non-continuing paths
			for l = 1:L
				if p.activePath[l]
					if keep[l,1] == false && keep[l,2] == false
						killPath(l, m, p)
					end
				end
			end
			
			# continue the other paths & duplicate if necessary
			for l = 1:L
				if keep[l,1] || keep[l,2]
					C = getArray_C(m+1, l, p)
					if keep[l,1] && keep[l,2]
						C[1, phi%2+1] = 0x00
						p.arrays_U[phi+1, l] = 0x00
						p.pathMetric[l, phi+1] = p.tmp_PM[l, 1]
						lprime = clonePath(l, m, phi, p)
						C = getArray_C(m+1, lprime, p)
						C[1, phi%2+1] = 0x01
						p.arrays_U[phi+1, lprime] = 0x01
						p.pathMetric[lprime, phi+1] = p.tmp_PM[l, 2]
					else # keep only one fork
						if keep[l,1]
							C[1, phi%2+1] = 0x00
							p.arrays_U[phi+1, l] = 0x00
							p.pathMetric[l, phi+1] = p.tmp_PM[l, 1]
						else
							C[1, phi%2+1] = 0x01
							p.arrays_U[phi+1, l] = 0x01
							p.pathMetric[l, phi+1] = p.tmp_PM[l, 2]
						end
					end
				end
			end
		end
		
		
		
		
		
	end
	
	function init_properties(p::Properties)
		# initialize data structures
		p.inactivePathIndices = collect(1:p.L)
		p.activePath = falses(p.L)

		
		@simd for lambda = 1:p.m+1
			@simd for s = 1:p.L
				p.arrayReferenceCount[lambda, s] = 0
			end
			p.inactiveArrayIndices[lambda] = collect(1:p.L)
		end
		
		@simd for i = 1:p.L
			@simd for j = 1:p.n
				p.pathMetric[i, j] = 0.0
			end
		end
	end
	
	
	
	# Binary Successive Cancelation List Decoder for Polar Codes (biAWGN channel)
	# input:
	#		p: Properties
	#		channel_LLRs: channel llrs ;)
	function scl_polar_decoder(p::Properties, channel_LLRs::Array{Float64,1})
		
		# local variables used quite often...
		n = p.n
		k = p.k
		m = p.m
		L = p.L
		frozen_bits = p.frozen_bits
		
		if p.listhandover == false
			init_properties(p) # initialize properties / re-set values
			
			## continue in line 2 of Alg. 12
			l = assignInitialPath(p)
			Ln = getArray_L(1, l, p)

			# init with probabilities from channel output
			@simd for b = 1:n
				Ln[b] = channel_LLRs[b];
			end
		else
			#rely on external initialization when listhandover == true
		end
		
		# main loop
		num_dbits = 0
		for phi = 0:n-1
			recursively_calc_L(m+1, phi, m, L, p)
			if haskey(frozen_bits, phi+1) # B_phi is a frozen bit
				@simd for l = 1:L
					if p.activePath[l]
						C = getArray_C(m+1, l, p)
						Ln = getArray_L(m+1, l, p)
						C[1, phi%2+1] = frozen_bits[phi+1]
						p.arrays_U[phi+1, l] = frozen_bits[phi+1]
						# update pathMetric
						if phi > 0
							p.pathMetric[l, phi+1] = calc_PM(p.pathMetric[l, phi], Ln[1], frozen_bits[phi+1])
						elseif p.listhandover
							p.pathMetric[l, phi+1] = calc_PM(p.pathMetric_handover[l], Ln[1], frozen_bits[phi+1])
						else
							p.pathMetric[l, phi+1] = 0.0
						end
					end
				end
			else
				continuePaths_UnfrozenBit(phi, m, L, p)
				num_dbits += 1
				if p.multipleCRC && num_dbits % (p.multipleCRC_interval + p.CRCbits) == 0
					# check CRC and prune list
					data = Array{UInt8,1}(undef, k+p.CRCbits)
					for l = 1:L
						if p.activePath[l]
							ind = 1
							for i = 1:phi+1
								if !haskey(frozen_bits, i)
									data[ind] = p.arrays_U[i,l]
									ind += 1;
								end
							end
							checksum = calc_CRC(data[1:ind-1-p.CRCbits], UInt(p.CRCbits), p.CRCpoly)
							if checksum != data[ind-p.CRCbits:ind-1] && sum(p.activePath) > 1
								# checksum does not match, remove path (but only if more then one path is left...)
								# TODO: this could be done better for FER reasons...
								killPath(l, m, p)
							end
						end
					end
				end
			end
			
			if phi % 2 == 1
				recursively_update_C(m+1, phi, m, L, p)
			end
		end
		
		# find the best codeword from the list
		if p.CRCbits > 0 && p.multipleCRC == false && p.listhandover == false # try to find codeword that fulfills CRC constraint
			t = copy(p.pathMetric[:,n])
			data = Array{UInt8,1}(undef, k+p.CRCbits)
			for l = 1:L
				lprime = findmin(t)[2]
				chat = getArray_C(1, lprime, p)[:,1]
				uhat = encode(chat)

				ind = 1
				for i = 1:n
					if !haskey(frozen_bits, i)
						data[ind] = uhat[i]
						ind += 1;
					end
				end
				
				#check CRC
				checksum = calc_CRC(data[1:k], UInt(p.CRCbits), p.CRCpoly)
				if checksum == data[k+1:end]
					# checksum matches :-)
					return chat
				end
				t[lprime] = Inf # set path metric to Inf, so it gets ignored in next iteration
			end
		end
		
		# no CRC match found, return first CW from list
		lprime = findmin(p.pathMetric[:,n])[2]
		C = getArray_C(1, lprime, p)[:,1]
		return C
	end	
end
