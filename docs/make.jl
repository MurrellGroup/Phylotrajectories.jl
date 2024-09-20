using Phylotrajectories
using Documenter

DocMeta.setdocmeta!(
    Phylotrajectories,
    :DocTestSetup,
    :(using Phylotrajectories);
    recursive = true,
)

makedocs(;
    modules = [Phylotrajectories],
    authors = "nossleinad <maximilian.danielsson@gmail.com> and contributors",
    sitename = "Phylotrajectories.jl",
    format = Documenter.HTML(;
        canonical = "https://MurrellGroup.github.io/Phylotrajectories.jl",
        edit_link = "main",
        assets = String[],
    ),
    pages = ["Home" => "index.md", "Simulate a count matrix" => "simulations.md"],
)

deploydocs(; repo = "github.com/MurrellGroup/Phylotrajectories.jl", devbranch = "main")
