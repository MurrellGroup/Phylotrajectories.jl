# `Phylotrajectories.jl` examples

This folder contains a self-contained OU-MCMC demo of the package.

## Files

| File                                | Description                                                                                                                                                                  |
|-------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `usage_example.ipynb`               | IJulia notebook walking through the full pipeline: load a wide-form count matrix, subsample, run OU-MCMC, build a HIPSTR consensus tree, plot diagnostic traces and a credibility-annotated tree. |
| `data/simulated_clone_data.csv`     | Wide-form simulated count matrix that backs the notebook — rows are clonotypes, columns are six cell-type subsets (`Subset1`…`Subset6`).                                       |

## Running the notebook

The notebook expects `IJulia` and a Julia ≥ 1.11 kernel:

```julia
using Pkg
Pkg.add("IJulia")
using IJulia; notebook(dir = "examples")
```

The first code cell calls `Pkg.activate(joinpath(@__DIR__, ".."))` so
that `using Phylotrajectories` resolves to the local checkout. If you've
already added the package via `Pkg.add(url = ...)` you can remove that
cell.

The notebook also reaches for `examples/data/simulated_clone_data.csv`
via `joinpath(@__DIR__, "data", ...)`, so it works out-of-the-box from
any clone of the repository.

## What the notebook produces

Outputs land in `examples/results/`:

```
OU_sim_<d>_<g>_<seed>_theta.pdf       # 6-panel diagnostic dashboard
OU_sim_<d>_<g>_hipster_<seed>_theta.pdf  # MCMC cloud + HIPSTR tree
HIPSTR_tree.newick                    # HIPSTR consensus (Newick)
HIPSTR_with_support.pdf               # HIPSTR with annotated supports
```
