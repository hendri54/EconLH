module ProductionFunctionsLH

using DocStringExtensions

export output, mproducts, validate_prod_fct, n_inputs, productivities, constant_returns
export CES, ces_output, ces_mproducts

abstract type AbstractProductionFunction{F1<: AbstractFloat} end

include("cobb_douglas.jl");
include("ces.jl");

end
