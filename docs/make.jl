using ActinFriction
using Documenter

DocMeta.setdocmeta!(ActinFriction, :DocTestSetup, :(using ActinFriction); recursive=true)

makedocs(;
    #modules=[ActinFriction],
    authors="Alexander Cumberworth <alex@cumberworth.org>",
    repo="https://github.com/cumberworth/ActinFriction.jl/blob/{commit}{path}#{line}",
    sitename="ActinFriction.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://cumberworth.github.io/ActinFriction.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/cumberworth/ActinFriction.jl",
    devbranch="master",
)
