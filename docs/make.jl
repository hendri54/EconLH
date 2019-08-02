using Documenter, EconLH

makedocs(
    modules = [EconLH],
    format = :html,
    checkdocs = :exports,
    sitename = "EconLH.jl",
    pages = Any["index.md"]
)

deploydocs(
    repo = "github.com/hendri54/EconLH.jl.git",
)
