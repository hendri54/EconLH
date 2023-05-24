module EconLH

# ToDo
# - Factor out production functions into a separate module

using DocStringExtensions
using Distributions, Random, StatsFuns

include("vectorLH.jl");
include("extreme_value_decision.jl");
include("production/productionFunctionsLH.jl");
include("gini.jl");
include("logistic.jl");


"""
    $(SIGNATURES)

Present value of a function of time.
First value is not discounted.
"""
function present_value(f, R :: T1, T :: T2) where
    {T1 <: AbstractFloat, T2 <: Integer}

    # return sum(f.(1 : T) ./ (R .^ (1:T))) .* R;
    pv = f(1);
    if T > one(T2)
        cumR = one(R);
        for t = 2 : T
            cumR *= R;
            pv += f(t) / cumR;
        end
    elseif T != one(T2)
        error("Invalid T: $T");
    end
    return pv
end


end # module
