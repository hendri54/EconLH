export ExtremeValueDecision, extreme_value_decision, extreme_value_decision_one, demeaned, pref_scale, draw_gumbel_shocks

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

"""
	$(SIGNATURES)

Use demeaned Gumbel preference shocks?
"""
demeaned(d :: ExtremeValueDecision{F}) where F = d.demeaned;

"""
	$(SIGNATURES)

Returns the scale parameter of the Gumbel preference shocks.
"""
pref_scale(d :: ExtremeValueDecision{F}) where F = d.prefScale;


## Decision: One alternative
function extreme_value_decision(d :: ExtremeValueDecision{F}, 
   value_iV :: AbstractVector{F}) where F

    return extreme_value_decision(value_iV, pref_scale(d); demeaned = demeaned(d));
end

function extreme_value_decision(value_iV :: AbstractVector{F}, prefScale :: F; 
   demeaned :: Bool = true) where F

    nTypes = length(value_iV);
    probV = ones(F, nTypes);
    valueV = copy(value_iV);
    if !demeaned
        valueV .+= (prefScale * F(EulerConst));
    end
    return probV, valueV
end


"""
   $(SIGNATURES)

Decision: Multiple alternatives.
Always returns a vector for `eVal_iV` (even with one agent).
   
OUT
   `prob_ixM`
      probability of choosing each option
   `eVal_iV`
      expected value for each type

TEST
   by simulation in `extreme_value_decision_test`
"""
function extreme_value_decision(d :: ExtremeValueDecision{F}, 
   value_ixM :: AbstractMatrix{F}) where  F <: Real

   return extreme_value_decision(value_ixM, pref_scale(d);
      demeaned = demeaned(d))
end

function extreme_value_decision(value_ixM :: AbstractMatrix{F}, prefScale :: F; 
   demeaned :: Bool = true) where F <: Real

   nTypes, nx = size(value_ixM);
   prob_ixM = similar(value_ixM);
   eVal_iV = Vector{F}(undef, nTypes);
   for j = 1 : nTypes
      prob_ixM[j,:], eVal_iV[j] = extreme_value_decision_one(vec(value_ixM[j,:]),
         prefScale; demeaned = demeaned);
   end
    return prob_ixM, eVal_iV
end

# Input is a Tuple with Vectors for each type (across all alternatives).
# Output is a Vector containing a Vector of decision probabilities for each type.
# test this +++++
function extreme_value_decision(value_ixV :: NTuple{nTypes, Vector{F}},
   prefScale :: F;  demeaned :: Bool = true) where {nTypes, F}

   # nx = length(first(value_ixV));
   prob_ixV = Vector{Vector{F}}(undef, nTypes);
   eVal_iV = Vector{F}(undef, nTypes);
   for j = 1 : nTypes
      prob_ixV[j], eVal_iV[j] = extreme_value_decision_one(value_ixV[j],
         prefScale; demeaned = demeaned);
   end
    return prob_ixV, eVal_iV

end


"""
	$(SIGNATURES)

Extreme value decision for one individual.

This is probably not efficient. See efficient implementation for `logsumexp` in `StatsFuns.jl`. But here we also need the scaled exp(valueV) directly for the probabilities.

# Arguments
- valueV: can be a Vector or a Tuple with a value for each alternative.
"""
function extreme_value_decision_one(valueV, prefScale :: F; 
   demeaned :: Bool = true) where F <: Real

   # Decision probability is log(sum(exp(V / prefScale)))
   # This needs to be nicely scaled to avoid overflow
   vMax = maximum(valueV) ./ prefScale .- F(4.0);

   # This uses a loop to ensure that probV is a Vector, even if valueV is a Tuple.
   # The following line is expensive
   probV = zeros(length(valueV));
   for (j, value) in enumerate(valueV)
      probV[j] = exp(value / prefScale - vMax);
   end

   # For each type: sum over alternatives
   expSum = sum(probV);

   # Prob of each choice
   probV ./= expSum;

   # Expected value
   eVal = prefScale * (vMax + log(expSum));

   if !demeaned
      eVal += prefScale * F(EulerConst);
   end

   return probV, eVal
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
function draw_gumbel_shocks(rng, prefScale :: F1, sizeV; demeaned :: Bool = true) where F1 <: Real

   drawM = prefScale .* rand(rng, Gumbel(zero(F1)),  sizeV...);
   if demeaned
       drawM .-= prefScale * F1(EconLH.EulerConst);
   end
   return drawM
end



# ----------------
