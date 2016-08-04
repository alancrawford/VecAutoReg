module VecAutoReg

using DataFrames, DataArrays, GLM, Distributions
import Base.show

typealias ArrayorDF  Union{Matrix{Float64},DataFrame}


include("varest.jl")

export 	readtable,	# Export this Data Frames Command to be able to load in data from REPL
		VarReg, 	# Constructor
		VarOutput, 	# VAR(p) Regression Output
		VarStable,  # Report whether the VAR(p) is stable 
		VarSoc,		# To directly call stats table
		VarOptLags, # Optimal Number of Lags Scored on AIC, BIC, HQIC
		VarFcast 	# Forecast 1 period ahead.
end