using Documenter
using ClusterAlgebras

DocMeta.setdocmeta!(ClusterAlgebras, :DocTestSetup, :(using ClusterAlgebras); recursive = true)

makedocs(;
    modules = [ClusterAlgebras],
    authors = "Benedikt Nagler <benedikt.nagler@protonmail.com>",
    sitename = "ClusterAlgebras.jl",
    format = Documenter.HTML(;
        canonical = "https://benedikt-nagler.github.io/ClusterAlgebras.jl",
        edit_link = "main",
        assets = String[],
        mathengine = Documenter.KaTeX(),
    ),
    pages = [
        "Home" => "index.md",
        "API" => "api.md",
    ],
    checkdocs = :exports,
)

deploydocs(; repo = "github.com/benedikt-nagler/ClusterAlgebras.jl", devbranch = "main")
