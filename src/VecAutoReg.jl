module VecAutoReg

using DataFrames, GLM, Distributions
import Base.show

include("varest.jl")

export 	readtable,	# Export this Data Frames Command to be able to load in data from REPL
		VarReg, 	# Constructor
		VarOutput, 	# VAR(p) Regression Output
		VarStable,  # Report whether the VAR(p) is stable 
		VarSoc,		# To directly call stats table
		VarOptLags, # Optimal Number of Lags Scored on AIC, BIC, HQIC
		VarFcast 	# Forecast 1 period ahead.
		GetSS_A		# Get State Space Matrix
end