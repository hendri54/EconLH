# EconLH

Economics related code.

# Discrete choice

This solves discrete choice problems subject to i.i.d. Gumbel preference shocks. An [`ExtremeValueDecision`](@ref) object can be constructed. Or [`extreme_value_decision`](@ref) can be called directly.

```@docs
ExtremeValueDecision
demeaned
pref_scale
extreme_value_decision
extreme_value_decision_one
draw_gumbel_shocks
```

# Production Functions

```@meta
CurrentModule = EconLH
```

```@docs
output
mproducts
n_inputs
productivities
constant_returns
validate_prod_fct
CES
ces_output
ces_mproducts
subst_elast
```

-------------