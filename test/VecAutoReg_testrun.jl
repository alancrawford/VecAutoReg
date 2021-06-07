#= Test for VAR code - using Stata var example =#
using VecAutoReg, CSV, DataFrames, StatsBase
cd(ENV["HOME"])
lutdata = CSV.read("./Git/VecAutoReg/testdata/lutkepohl12.csv", DataFrame);
constant_term = false; # Set Constant term

#= ************** Call VAR from Array ************** =#
Tobs, NumCols = size(lutdata);
K = 3;
X = zeros(Tobs-1,K);
X[:,1] = lutdata[2:end, :dln_inv];
X[:,2] = lutdata[2:end, :dln_inc];
X[:,3] = lutdata[2:end, :dln_consump];

Pmax = 10;
VarOptLags(X, Pmax, constant_term)
var2 = VarReg(X, 2, constant_term);
VarOutput(var2)
VarStable(var2)

#= ****************** VAR Data Frame Call ******************* =#
TSdata = lutdata[2:end,:];
eqvars = names(lutdata)[6:2:10];

Pmax = 10;
VarOptLags(lutdata, Pmax, eqvars, constant_term )
var2 = VarReg(TSdata, 2,  eqvars, constant_term);
VarOutput(var2)
VarStable(var2)




