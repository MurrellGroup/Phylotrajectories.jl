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

    warnonly = [:missing_docs],
    pages = [
        "Home" => "index.md",
        "Getting started" => "getting_started.md",
        "Model & parameters" => "models.md",
        "Simulating a count matrix" => "simulations.md",
        "Worked example" => "example_workflow.md",
        "API reference" => "api.md",
    ],
)

deploydocs(; repo = "github.com/MurrellGroup/Phylotrajectories.jl", devbranch = "main")
