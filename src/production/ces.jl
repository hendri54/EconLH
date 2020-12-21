"""
  CES

CES aggregator of the form

`Y = A * [sum of alpha (X ^ rho)] ^ (1/rho)`
`Q = sum of (alpha X) ^ rho`

Sum of `alpha` = 1. Which is why the `alpha` are not taken to the power `rho`.

Accommodates matrix inputs (T x N)

Best to store AV in the object
- little cost in terms of efficiency
- consistency ensured
"""
mutable struct CES{F1} <: AbstractProductionFunction{F1}
  substElast :: F1
  # relative weights; Nx1
  alphaV :: Vector{F1}
  # Productivities
  AV :: Vector{F1}
end

"""
	$(SIGNATURES)

Substitution elasticity.
"""
subst_elast(fS :: CES{F1}) where F1 = fS.substElast;

"""
	$(SIGNATURES)

Factor neutral productivities.
"""
productivities(fS :: CES{F1}) where F1 = fS.AV;

"""
	$(SIGNATURES)

Number of inputs.
"""
n_inputs(fS :: CES{F1}) where F1 = length(fS.alphaV);

"""
	$(SIGNATURES)

Does the production function have constant returns to scale?
"""
constant_returns(fS :: CES{F1}) where F1 = true;

# function CES(sElast :: F1, aV :: Vector{F1}, prodV :: Vector{F1}) where F1 <: AbstractFloat

#     fS = new();
#     set_params(fS, substElast = sElast, alphaV = aV, AV = prodV);
#     return fS
#   end
# end


"""
	$(SIGNATURES)

Validate a CES production function.
"""
function validate_prod_fct(this :: CES{F1}) where F1 <: AbstractFloat
    # T, n = size(fS.xM);
    isValid = check_subst_elast(this.substElast)
    # check more +++
end

function check_subst_elast(substElast :: F1) where 
  F1 <: AbstractFloat
    return (substElast > 0.01  &&  substElast < 50)
end


## Set parameters in CES
# function set_params(fS :: CES;  substElast = [], alphaV = [], AV = [])
#   # Assign, then check consistency
#   if !isempty(substElast)
#     fS.substElast = substElast;
#   end
#   if !isempty(alphaV)
#     fS.alphaV = alphaV;
#   end
#   if !isempty(AV)
#     fS.AV = AV;
#   end
#   # if !isempty(xM)
#   #   fS.xM = xM;
#   # end

#   # Derived
#   fS.rho = 1.0 - 1.0 / fS.substElast;

#   check(fS);
# end


# # ----  return rho from subst elast
function curvature(fS :: CES{F1}) where F1 <: AbstractFloat
  return one(F1) / fS.substElast - one(F1);
end


## -------  Returns sum(alpha * x ^ rho) ^ (1/rho)
# With rho = 0 (Cobb-Douglas) this returns 1.0.
function ces_q(alphaV :: Vector{F1}, xM :: Matrix{F1}, rho :: F1) where F1 <: AbstractFloat
  T = size(xM, 1);
  qV = zeros(F1, T);
  for t = 1 : T
      # Careful here: alphaV * xM could become a matrix
      qV[t] = sum(alphaV .* (xM[t,:] .^ rho));
  end
  return qV
end


"""
	$(SIGNATURES)

Output given inputs. If substitution elasticity is close to 1, switch to Cobb-Douglas.
"""
function output(fS :: CES{F1}, xM :: Matrix{F1}) where F1 <: AbstractFloat
  qV = ces_q(fS.alphaV, xM, curvature(fS));
  # This does not overflow when curvature = 1.0
  yV = productivities(fS) .* (qV .^ (one(F1) / curvature(fS)));
  if switch_cobb_douglas(subst_elast(fS), yV)
    # Handle numerical overflow by switching to Cobb-Douglas
    yV = productivities(fS) .* output_cobb_douglas(xM, fS.alphaV);
  end
  return yV
end

# Decide whether Cobb-Douglas needs to be used
function switch_cobb_douglas(substElast, xV)
  if abs(substElast - 1.0) < 0.01
    return true
  elseif (abs(substElast - 1.0) < 0.05)  &&  any(isinf.(xV))
    return true
  else
    return false
  end
end


"""
	$(SIGNATURES)

Marginal products (TxN).
"""
function mproducts(fS :: CES{F1}, xM :: Matrix{F1}) where F1 <: AbstractFloat

  T = length(fS.AV)
  N = length(fS.alphaV);
  rho = curvature(fS);
  # qV is T x 1
  qV = ces_q(fS.alphaV, xM, rho);
  aqV = fS.AV .* (qV .^ (one(F1) / rho - one(F1)));
  # mpM = similar(xM);
  # for j = 1 : N
  #   mpM[:,j] = aqV .* (fS.alphaV[j] .^ rho) .* (xM[:,j] .^ (rho - 1));
  # end
  if switch_cobb_douglas(subst_elast(fS), aqV)
    mpM = fS.AV .* mproducts_cobb_douglas(xM, fS.alphaV);
  else
    mpM = aqV .* (xM .^ (rho-one(F1))) .* fS.alphaV';
  end
  # @assert (all(abs.(mp2M - mpM) .< 1e-5))
  return mpM
end


function make_test_ces(T, n, substElast)
  alphaV = collect(range(1.0, 2.0, length = n));
  alphaV = alphaV ./ sum(alphaV);
  AV = collect(LinRange(0.5, 9.0, T));

  fS = CES(substElast, alphaV, AV);
  @assert validate_prod_fct(fS)
  return fS
end

# ------------