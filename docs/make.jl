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
        "Tutorial" => "tutorial.md",
        "Manual" => [
            "Quivers" => "quivers.md",
            "Seeds and mutation" => "seeds.md",
            "Coefficients" => "coefficients.md",
            "Types and classification" => "classification.md",
            "Mutation classes" => "classes.md",
            "Green sequences and DT" => "green.md",
            "Friezes" => "friezes.md",
            "Grassmannians" => "grassmannian.md",
            "Interchange formats" => "interop.md",
        ],
        "Errors" => "errors.md",
        "API index" => "api.md",
    ],
    checkdocs = :exports,
)

deploydocs(; repo = "github.com/benedikt-nagler/ClusterAlgebras.jl", devbranch = "main")
