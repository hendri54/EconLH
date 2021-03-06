module EconLH

# ToDo
# - Factor out production functions into a separate module

using DocStringExtensions
using Distributions, Random, StatsFuns

include("vectorLH.jl");
include("extreme_value_decision.jl");
include("production/productionFunctionsLH.jl");


"""
    $(SIGNATURES)

Present value of a function of time.
First value is not discounted.
"""
function present_value(f, R :: T1, T :: T2) where
    {T1 <: AbstractFloat, T2 <: Integer}

    return sum(f.(1 : T) ./ (R .^ (1:T))) .* R;
end


end # module
