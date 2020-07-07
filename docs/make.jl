using Documenter, EconLH, EconLH.LatexLH

makedocs(
    modules = [EconLH],
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    checkdocs = :exports,
    sitename = "EconLH.jl",
    pages = Any["index.md"]
)

# ---------