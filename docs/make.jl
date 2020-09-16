Pkg.activate("./docs");

using Documenter, EconLH, FilesLH

makedocs(
    modules = [EconLH],
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    checkdocs = :exports,
    sitename = "EconLH.jl",
    pages = Any["index.md"]
)

pkgDir = rstrip(normpath(@__DIR__, ".."), '/');
@assert endswith(pkgDir, "EconLH")
deploy_docs(pkgDir);

Pkg.activate(".");

# ---------