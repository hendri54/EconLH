# EconLH

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

# LatexLH module

```@meta
CurrentModule = EconLH.LatexLH
```

This contains Latex related code. 

## Writing Beamer Slides

```@docs
write_figure_slides
figure_slide
```

## Writing Latex Tables

```@docs
CellColor
Cell
Table
color_string
cell_string
nrows
add_row!
make_row
write_table
```

-------------