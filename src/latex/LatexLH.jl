"""
Latex code

For notation: set up a `Dict{String, String}`. No need for code.
"""
module LatexLH

using DocStringExtensions

export write_figure_slide, figure_slide

include("table.jl");
include("parameter_table.jl")
include("beamer.jl")
# include("latex/symbol_table.jl")

end
