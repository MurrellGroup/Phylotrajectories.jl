```@meta
CurrentModule = Phylotrajectories
```

# Phylotrajectories.jl

[Phylotrajectories.jl](https://github.com/MurrellGroup/Phylotrajectories.jl)
infers **cell-type phylogenies** from single-cell clonotype-by-cell-type
count matrices. It treats clonotype frequencies as an **Ornstein–Uhlenbeck
(OU) process** — a Brownian motion with mean-reversion towards an
equilibrium — diffusing along an unknown tree of cell phenotypes, and
samples the joint posterior over topology, branch lengths and OU
parameters via Metropolis-Hastings MCMC. See [`OUContinuousModel`](@ref)
for the full specification.

The package also provides:

- [`sim_count_matrix`](@ref) for **simulating** count matrices under a
  diffusion process on a tree, used for validation and for building
  tutorials.
- [`run_ou_and_build_clone_matrix`](@ref) to convert a posterior tree
  ensemble into a per-clonotype × per-node frequency matrix.
- A small set of plotting helpers (`PlotTreeOnUmap`,
  `PlotTreeOnUmapNoAnim`, `PlotTreeOnUmapNoAnimShadow`) for projecting
  inferred trees onto a 2-D embedding such as a UMAP of the cells.
- [`recombine`](@ref) for merging two cell types in a count matrix.

## Installation

```julia
using Pkg
Pkg.add(url = "https://github.com/MurrellGroup/Phylotrajectories.jl")
```

## A 30-second tour

```julia
using Phylotrajectories

# 1. read a long-form per-cell table (Barcode, cell_types, Clonotype, …)
clono_info, cluster_names, cluster_sizes, count_matrix = import_count_matrix(
    "examples/data/simulated_clone_data.tsv",
    :Clonotype, :cell_types, :TRB_cdr3aa,
)

# 2. fit the OU model with MCMC
plot_init, init_tree, trees, LLs, models, root_ps, upd =
    tree_inference(
        OUContinuousModel(n_samples = 200, burn_in = 2_000, sample_interval = 50),
        cluster_names, count_matrix;
        eqmu = 1.5, eqtheta = 0.1, v = 1.0, d = 0.5, g = 0.5,
    )
```

## Documentation pages

- **[Getting started](getting_started.md)** — installation, importing data, running inference.
- **[Model & parameters](models.md)** — full reference for `OUContinuousModel` and its sampler primitives.
- **[Simulating data](simulations.md)** — how to use `sim_count_matrix` to validate a workflow.
- **[Worked example](example_workflow.md)** — end-to-end OU-MCMC + UMAP overlay on simulated data, mirroring the [`examples/usage_example.ipynb`](https://github.com/MurrellGroup/Phylotrajectories.jl/blob/main/examples/usage_example.ipynb) notebook.
- **[API reference](api.md)** — every exported symbol with its docstring.

## Citation

If you use Phylotrajectories.jl in academic work, please cite the
[GitHub repository](https://github.com/MurrellGroup/Phylotrajectories.jl).
