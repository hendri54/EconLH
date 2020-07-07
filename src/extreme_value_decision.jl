export ExtremeValueDecision, extreme_value_decision, demeaned, pref_scale, draw_gumbel_shocks

const EulerConst = 0.5772;

"""
   $(SIGNATURES)

Solve a decision problem with type I extreme value shocks

# Arguments
- value_ixM
      values of alternatives
      rows are types; columns are alternatives
      works with a single type
- prefScale
      scale of type I extreme value shocks
"""
struct ExtremeValueDecision{F <: AbstractFloat}
    prefScale :: F
    # Preference shocks demeaned?
    demeaned :: Bool
    dbg :: Bool
end

function ExtremeValueDecision(prefScale :: F, dbg :: Bool = false) where
    F <: AbstractFloat

    return ExtremeValueDecision(prefScale, true, dbg)
end

demeaned(d :: ExtremeValueDecision{F}) where F = d.demeaned;
pref_scale(d :: ExtremeValueDecision{F}) where F = d.prefScale;


## Decision: One alternative
function extreme_value_decision(d :: ExtremeValueDecision{F}, value_iV :: Vector{F}) where
    F <: AbstractFloat

    nTypes = length(value_iV);
    probV = ones(F, nTypes);
    valueV = copy(value_iV);
    if !demeaned(d)
        valueV .+= (pref_scale(d) * F(EulerConst));
    end
    return probV, valueV
end


"""
   $(SIGNATURES)

Decision: Multiple alternatives.
Always returns a vector for `eVal_iV` (even with one agent).
   
OUT
   prob_ixM
      probability of choosing each option
   eVal_iV
      expected value for each type

TEST
   by simulation in extreme_value_decision_test
"""
function extreme_value_decision(d :: ExtremeValueDecision{F}, 
   value_ixM :: Matrix{F}) where  F <: AbstractFloat

    nTypes, nx = size(value_ixM);

    # Decision probability is log(sum(exp(V / prefScale)))
    # This needs to be nicely scaled to avoid overflow
    vMax_iV = maximum(value_ixM ./ pref_scale(d), dims = 2) .- 4.0;
    # The following line is expensive
    exp_ixM = exp.(value_ixM ./ pref_scale(d) .- vMax_iV);

    # For each type: sum over alternatives
    expSum_iV = sum(exp_ixM, dims = 2);

    # Prob of each choice
    prob_ixM = exp_ixM ./ expSum_iV;

    # Expected value
    eVal_iV = vec(pref_scale(d) .* (vMax_iV + log.(expSum_iV)));

    if !demeaned(d)
      eVal_iV .+= pref_scale(d) * EulerConst;
    end

    if d.dbg
       @assert isreal(exp_ixM)
       @assert size(exp_ixM) == size(value_ixM)
       @assert size(prob_ixM) == size(value_ixM)
       @assert size(eVal_iV) == (nTypes,)
    end

    return prob_ixM, eVal_iV
end

# """
# Same with vector input (b/c a vector is not the same as an array with one row)
# """
# function extreme_value_decision(valueV :: Array{T,1},  prefScale :: T,
#     dbg :: Bool)   where T <: Real
#
#     return extreme_value_decision(collect(valueV'), prefScale, dbg)
# end


"""
	$(SIGNATURES)

Draw Gumbel shocks, optionally demeaned. Returns vector, even if `sizeV == 1`.
"""
function draw_gumbel_shocks(rng, prefScale :: F1, sizeV; demeaned :: Bool = true) where F1 <: AbstractFloat

   drawM = prefScale .* rand(rng, Gumbel(zero(F1)),  sizeV...);
   if demeaned
       drawM .-= prefScale * F1(EconLH.EulerConst);
   end
   return drawM
end



# ----------------
