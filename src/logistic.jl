## -------------  Logistic function

export GeneralizedLogistic, logistic;

"""
	$(SIGNATURES)

Generalized logistic function object.

``
lb + (ub - lb) / (1 + f0 * exp(-slope * x)) ^ twister
``

Note that slope can be set to 1 if the input is `x / x0`. 
If two logistic functions have the same `f0` and bounds, they have the same `f(0)`.

# Parameters
- `lb` and `ub` are the bounds of `f(x)`.
- `f0 >= 0` roughly shifts the curve left/right. Higher `f0` implies lower f(0) or shifts the curve right.
- `slope > 0` is the slope.
- `twister > 0` produces asymmetry (`twister = 1` is symmetric).
"""
Base.@kwdef mutable struct GeneralizedLogistic{F1}
    lb :: F1 = zero(F1)
    ub :: F1 = one(F1)
    f0 :: F1 = one(F1)
    slope :: F1 = one(F1)
    twister :: F1 = one(F1)
end


"""
	$(SIGNATURES)

Generalized logistic function.
"""
function logistic(x; lb = 0.0, ub = 1.0, f0 = 1.0, slope = 1.0, twister = 1.0)
    return lb .+ (ub - lb) ./ ((1 .+ f0 .* exp.(-slope .* x)) .^ twister);
end

function logistic(glf :: GeneralizedLogistic{F1}, x) where F1
    return logistic(x; lb = glf.lb, ub = glf.ub, f0 = glf.f0,
        slope = glf.slope, twister = glf.twister)
end

# -----------