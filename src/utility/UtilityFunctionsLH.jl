module UtilityFunctionsLH

abstract type UtilityFunction end
abstract type UtilityOneArg <: UtilityFunction end

include("crra.jl")

end
