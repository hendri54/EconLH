module ProductionFunctionsLH

using DocStringExtensions

export validate_prod_fct, n_inputs, productivities, constant_returns
export CES, output, mproducts

abstract type AbstractProductionFunction{F1<: AbstractFloat} end

include("cobb_douglas.jl");
include("ces.jl");

end
