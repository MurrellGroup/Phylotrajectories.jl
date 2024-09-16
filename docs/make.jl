using Phylotrajectories
using Documenter

DocMeta.setdocmeta!(Phylotrajectories, :DocTestSetup, :(using Phylotrajectories); recursive=true)

makedocs(;
    modules=[Phylotrajectories],
    authors="nossleinad <maximilian.danielsson@gmail.com> and contributors",
    sitename="Phylotrajectories.jl",
    format=Documenter.HTML(;
        canonical="https://nossleinad.github.io/Phylotrajectories.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/nossleinad/Phylotrajectories.jl",
    devbranch="main",
)
