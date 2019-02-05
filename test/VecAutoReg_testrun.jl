#= Test for VAR code - using Stata var example =#
using VecAutoReg, CSVFiles, DataFrames
lutdata = DataFrame(load("./Git/VecAutoReg/testdata/lutkepohl12.csv"));
constant_term = true; # Set Constant term

#= ************** Call VAR from Array ************** =#
Tobs, NumCols = size(lutdata);
K = 3;
X = zeros(Tobs-1,K);
X[:,1] = lutdata[:dln_inv][2:end,:];
X[:,2] = lutdata[:dln_inc][2:end,:];
X[:,3] = lutdata[:dln_consump][2:end,:];

Pmax = 10;
VarOptLags(X, Pmax, constant_term)
var2 = VarReg(X, 2, constant_term);
VarOutput(var2)
VarStable(var2)

#= ****************** VAR Data Frame Call ******************* =#
TSdata = lutdata[2:end,:];
eqvars = names(lutdata)[6:2:10];

Pmax = 10;
VarOptLags(TSdata, Pmax, eqvars, constant_term )
var2 = VarReg(TSdata, 2,  eqvars, constant_term);
VarOutput(var2)
VarStable(var2)




