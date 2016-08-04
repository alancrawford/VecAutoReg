#= VAR code =#

type VarReg
	T			:: Int64				# Number of time series observations
	N 			:: Int64 				# Number of series
	lags 		:: Int64 				# Number of lags
	X 			:: ArrayorDF 			# Data 
	eqnames 	:: Vector{Symbol}		# Variable names
	consterm	:: Bool					# Boolean for inclusion of a constant term
	Bhat 		:: Matrix{Float64}		# Estimates of parameters
	SigmaU		:: Matrix{Float64} 		# Covariance Matrix 
	CovBhat		:: Matrix{Float64} 		# Covariance Matrix of Bhat
	se			:: Matrix{Float64}		# Standard Errors	
	err_dist	:: Distributions.MvNormal 	# Fitted dist of residuals

	function VarReg(X::Matrix{Float64}, p::Int64,consterm::Bool=true; varname::ASCIIString="x")
		#= Setup =#
		(T,N) = size(X)
		is(consterm,true) ? 
			@assert(T*N > p*N*N, "Not enough data points to estimate $p lags") :
			@assert(T*N > p*N*(N+1), "Not enough data points to estimate $p lags")
		Y = X[p+1:T,:]
		Z = hcat(ones(T-p), zeros(T-p,N*p))
		for q = 1:p
			Z[:,N*(q-1)+2:N*q+1] = X[p-q+1:end-q,:] 
		end
		is(consterm,true) ? nothing :	Z = Z[:,2:end]  # Omit column of ones if constermtant is false

		#= Estimation =#
		Ttilde = T-p
		invZZ = inv(Z'*Z)
		Bhat = Z\Y
		Uhat = Y - Z*Bhat
		SigmaU = Uhat'*Uhat/Ttilde
		Gamma = Z'*Z/Ttilde
		CovBhat = kron(invZZ,SigmaU)
		is(consterm,true) ?
			se = reshape(sqrt(diag(CovBhat)),N,N*p+1)' :
			se = reshape(sqrt(diag(CovBhat)),N,N*p)'

		# Because no Data Frame create series names
		eqnames = [symbol(varname*"$i") for i in 1:N]
		err_dist = MvNormal(zeros(N),SigmaU)
	
		new(T, N, p, X, eqnames , consterm, Bhat, SigmaU, CovBhat, se, err_dist)
	end

	function VarReg(df::DataFrame, p::Int64,  eqvars::Vector{Symbol}, consterm::Bool=true)
		#= Setup =#
		X = df[:,eqvars]
		(T,N) = size(X)
		is(consterm,true) ? 
			@assert(T*N > p*N*N, "Not enough data points to estimate $p lags") :
			@assert(T*N > p*N*(N+1), "Not enough data points to estimate $p lags")
		ctr = 1
		Y = Array(Float64,T-p,N)
		for i in 1:N
			Y[:,ctr] = X[eqvars[ctr]][p+1:end]
			ctr +=1
		end
		Z = [ones(T-p) zeros(T-p,N*p)]
		ctr = 1
		for q = 1:p
			for i in 1:length(eqvars)
				ctr += 1
				Z[:,ctr] = X[eqvars[i]][p-q+1:end-q] 
			end
		end
		is(consterm,true) ? nothing :	Z = Z[:,2:end]  # Omit column of ones if constant is false

		#= Estimation =#
		Ttilde = T-p
		invZZ = inv(Z'*Z)
		Bhat = Z\Y
		Uhat = Y - Z*Bhat
		SigmaU = Uhat'*Uhat/Ttilde
		Gamma = Z'*Z/Ttilde
		CovBhat = kron(invZZ,SigmaU)
		is(consterm,true) ?
			se = reshape(sqrt(diag(CovBhat)),N,N*p+1)' :
			se = reshape(sqrt(diag(CovBhat)),N,N*p)'
	
		new(T, N, p, X, eqvars , consterm, Bhat, SigmaU, CovBhat, se)
	end

end

function show(io::IO, vrp::VarReg)
	msg = "Ran VAR with:\n"
	msg *= "\t- $(vrp.T) time periods\n "
	msg *= "\t- $(vrp.N) series\n "
	msg *= "\t- $(vrp.lags) lags\n "
	msg *= "\t- Constant Term: $(vrp.consterm)\n"
	msg *= "\t- NumObs used in estimation = $((vrp.T-vrp.lags)*vrp.N)\n"
	print(io,msg)
end

function VarOutput( vrp :: VarReg)
	pK2 = vrp.lags*vrp.N*vrp.N
	model_df = vrp.T*vrp.N-pK2-vrp.N
	zstat = vrp.Bhat./vrp.se
	pval = ccdf(FDist(1,model_df),abs2(zstat))
	ci_low = vrp.Bhat - 1.96vrp.se
	ci_high = vrp.Bhat + 1.96vrp.se

	# Add pretty output
	is(vrp.consterm,true) ? RowNames = ["_Cons"] : RowNames = []
	for q = 1:vrp.lags
		for w = 1:vrp.N
			push!(RowNames,"L$q.$(vrp.eqnames[w])")
		end
	end
	ColNames = ["Coef. ","  SE  ","  z  "," P>|z| ","[CI 5%,","CI 95%]"]
	for w in 1:vrp.N
		println("\n ** Eq: $(vrp.eqnames[w]) **\n")
		regout = round([vrp.Bhat[:,w] vrp.se[:,w] zstat[:,w] pval[:,w] ci_low[:,w] ci_high[:,w]],4)
		show(CoefTable(regout, ColNames,  RowNames))
	end
end

function VarStable(vrp::VarReg)
	A = [vrp.Bhat[2:end,:]'; [eye(vrp.N*(vrp.lags-1)) zeros(vrp.N*(vrp.lags-1),vrp.N)]]
	MaxModEigVal = abs(eigvals(A))[1];
	EigOut = round(MaxModEigVal,4)
	<(MaxModEigVal,1.0) ?
		println("Max Eigenvalue is $EigOut < 1. VAR stabiliity condition is met." ) :
		println("Max Eigenvalue is $EigOut > 1. VAR stabiliity condition is not met.")
end

function VarSoc(vrp::VarReg)
	pK2 = vrp.lags*vrp.N*vrp.N
	Ttilde = vrp.T - vrp.lags
	LL = log(det(vrp.SigmaU))
	AIC =  LL + 2*pK2/Ttilde
	SBIC = LL + pK2*log(Ttilde)/Ttilde
	HQIC = LL + 2*pK2*log(log(Ttilde))/Ttilde
	return [AIC SBIC HQIC]
end

function VarOptLags( X::Matrix{Float64} , Pmax::Int64, consterm::Bool=true)
	soc = Array(Float64,Pmax,3)
	for q in 1:Pmax
		vrp = VarReg(X, q, consterm)
		soc[q,:] = VarSoc(vrp)
	end
	MinCrit = zeros(Float64,3)
	for j in 1:3
		 (),MinCrit[j]  = findmin(soc[:,j])
	end
	ColNames = ["  AIC ","  SBIC "," HQIC "];
	RowNames  = ["lag = $i" for i in 1:Pmax];
	CritRowNames = ["Opt. Lags"]
	println("-----------------------------------\n")
	println("\t** Optimal Lag Length **\n")
	println("-----------------------------------\n")
	show(CoefTable(MinCrit', ColNames, CritRowNames))
	println("-----------------------------------\n")	
	show(CoefTable(soc, ColNames, RowNames))
	return MinCrit'
end

function VarOptLags(X::DataFrame, Pmax::Int64,  eqvars::Vector{Symbol}, consterm::Bool=true)
	soc = Array(Float64,Pmax,3)
	for q in 1:Pmax
		vrp = VarReg(X, q, consterm, eqvars)
		soc[q,:] = VarSoc(vrp)
	end
	MinCrit = zeros(Float64,3)
	for j in 1:3
		 (),MinCrit[j]  = findmin(soc[:,j])
	end
	ColNames = ["  AIC ","  SBIC "," HQIC "];
	RowNames  = ["lag = $i" for i in 1:Pmax];
	CritRowNames = ["Opt. Lags"]
	println("-----------------------------------\n")
	println("\t** Optimal Lag Length **\n")
	println("-----------------------------------\n")
	show(CoefTable(MinCrit', ColNames, CritRowNames))
	println("-----------------------------------\n")	
	show(CoefTable(soc, ColNames, RowNames))
	return MinCrit'
end


function VarFcast(vrp::VarReg, X::Array{Float64,1})
	T, N = size(X)
	@assert(==(N,vrp.N), "Input Data not conformable. Should be T x $vrp.N")
	==(vrp.consterm,true) ? 
		Xprime = [ones(Float64, T) X]*vrp.Bhat : 
		Xprime =  X*vrp.Bhat	
	return Xprime
end

function VarFcast(vrp::VarReg, X::Array{Float64,2})
	T, N = size(X)
	@assert(==(N,vrp.N), "Input Data not conformable. Should be T x $vrp.N")
	==(vrp.consterm,true) ? 
		Xprime = [ones(Float64,T) X]*vrp.Bhat : 
		Xprime =  X*vrp.Bhat	
	return Xprime
end



