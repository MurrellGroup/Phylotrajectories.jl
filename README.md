# Phylotrajectories

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://MurrellGroup.github.io/Phylotrajectories.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://MurrellGroup.github.io/Phylotrajectories.jl/dev/)
[![Build Status](https://github.com/MurrellGroup/Phylotrajectories.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/MurrellGroup/Phylotrajectories.jl/actions/workflows/CI.yml?query=branch%3Amain)

`Phylotrajectories.jl` infers cell-type phylogenies from single-cell
clonotype-by-cell-type count matrices. It treats clonotype frequencies
as an **Ornstein–Uhlenbeck (OU) process** — a Brownian motion with
mean-reversion towards an equilibrium — diffusing along an unknown tree
of cell phenotypes, and samples the joint posterior over topology,
branch lengths and OU parameters via Metropolis-Hastings MCMC.

See the [docs](https://MurrellGroup.github.io/Phylotrajectories.jl/dev/) for
the full reference.

## Installation

```julia
using Pkg
Pkg.add(url = "https://github.com/MurrellGroup/Phylotrajectories.jl")
```

## Quick start on the bundled simulated dataset

The repository ships with a tiny simulated dataset under `examples/data/`
and a runnable end-to-end notebook at
[`examples/usage_example.ipynb`](examples/usage_example.ipynb).

```julia
using Phylotrajectories

clono_info, cluster_names, cluster_sizes, count_matrix = import_count_matrix(
    "examples/data/simulated_clone_data.tsv",
    :Clonotype, :cell_types, :TRB_cdr3aa,
)

plot_init, init_tree, trees, LLs, models, root_ps, upd =
    tree_inference(
        OUContinuousModel(burn_in = 2_000, sample_interval = 50, n_samples = 200),
        cluster_names, count_matrix;
        eqmu = 1.5, eqtheta = 0.1, v = 1.0, d = 0.5, g = 0.5,
    )

ladderize!.(trees)
hip, node2logcred, node2support = HIPSTR(trees; getcred = true, getsupport = true)
```

The notebook adds the HIPSTR credibility plot, the per-clone × per-node
frequency matrix produced by `run_ou_and_build_clone_matrix`, and a UMAP
overlay built with `PlotTreeOnUmap*`.

## Importing data

`import_count_matrix` accepts two file shapes:

### Long-form (one row per cell)

```julia
clono_info, cluster_names, cluster_sizes, count_matrix = import_count_matrix(
    "data/clone_data_HDM.tsv",
    :Clonotype, :cell_types, :TRB_cdr3aa,
    cluster_filters = ["Proliferating"],
)
```

`:Clonotype`, `:cell_types` and `:TRB_cdr3aa` are the column names that hold
the clonotype, cell-type label and TRB CDR3 sequence respectively. The
`cluster_filters` keyword drops named clusters before pivoting.

### Wide-form (already-pivoted CSV)

```julia
clono_info, cluster_names, cluster_sizes, count_matrix =
    import_count_matrix("data/Clone_counts_HDM.csv")
```

## Performing inference

```julia
plot_init, init_tree, trees, LLs, models, root_ps, upd =
    tree_inference(
        OUContinuousModel(burn_in = 2_000, sample_interval = 50, n_samples = 200),
        cluster_names, count_matrix;
        eqmu = 1.5, eqtheta = 0.1, v = 1.0, d = 0.5, g = 0.5,
    )
```

The returned tuple is:

- `plot_init`  — `Plots.Plot` of the starting tree (handy for sanity checks),
- `init_tree`  — `FelNode` actually used to start MCMC,
- `trees`      — `Vector{FelNode}` of posterior topology samples,
- `LLs`        — log-likelihoods at the sampled points,
- `models`     — sampled `(θ, v, μ)` triples,
- `root_ps`    — sampled root state distributions,
- `upd`        — the `Update` object holding acceptance ratios and proposal stats.

See the [Models & parameters](https://MurrellGroup.github.io/Phylotrajectories.jl/dev/models/)
reference for every available knob on `OUContinuousModel`.

## Repository layout

```
src/                # package source
  inference/        # OU MCMC samplers
  viz/              # tree-on-UMAP plotting helpers
  utils/            # post-processing helpers (HIPSTR, clone matrices, metrics)
  importing.jl      # import_count_matrix
  simulations.jl    # sim_count_matrix
docs/               # Documenter source
examples/           # bundled simulated dataset + Jupyter notebook
test/               # unit tests
```
