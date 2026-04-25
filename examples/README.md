# `Phylotrajectories.jl` examples

This folder contains a self-contained demo of the package on a small
simulated dataset.

## Files

| File                                  | Description                                                                                                                                                  |
|---------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `data/simulated_clone_data.tsv`       | Long-form per-cell table (Barcode, `cell_types`, `Clonotype`, CDR3 columns). 4&nbsp;000 cells across 8 cell types and 40 clonotypes.                          |
| `data/simulated_UMAP_coords.tsv`      | Per-cell 2-D UMAP-style embedding (Barcode, `cell_types`, `UMAP_1`, `UMAP_2`).                                                                                |
| `generate_simulated_data.jl`          | Julia script that produces both files with `Phylotrajectories.sim_count_matrix` (Brownian-motion-on-a-tree count simulation, plus a per-clonotype "bias"). |
| `usage_example.ipynb`                 | IJulia notebook walking through the full OU-MCMC + HIPSTR + UMAP pipeline on the simulated data.                                                              |

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

## Regenerating the simulated data

```bash
julia --project=. examples/generate_simulated_data.jl
```

The script is deterministic — `Random.seed!(20260425)` at the top — so
re-running it reproduces the bundled TSVs exactly.
