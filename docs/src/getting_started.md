```@meta
CurrentModule = Phylotrajectories
```

# Getting started

## Installation

The package isn't yet registered in the Julia General registry. Install
it directly from the GitHub repository:

```julia
using Pkg
Pkg.add(url = "https://github.com/MurrellGroup/Phylotrajectories.jl")
```

## Importing data

Phylotrajectories works on a `cell_types × clonotypes` integer count
matrix. It can be built from two file shapes:

### Long-form (one row per cell)

A TSV with at least three columns: a clonotype label, a cell-type label,
and a TRB CDR3 amino-acid sequence (extra columns are tolerated).

```julia
using Phylotrajectories

clono_info, cluster_names, cluster_sizes, count_matrix = import_count_matrix(
    "data/clone_data_HDM.tsv",
    :Clonotype,        # column with clonotype IDs
    :cell_types,       # column with cell-type labels
    :TRB_cdr3aa,       # column with TRB CDR3 sequences
    cluster_filters = ["Proliferating"],   # cell types to ignore
)
```

`clono_info` is a `DataFrame` mapping each cell type to the set of TRB
CDR3 strings observed in it. `cluster_sizes` maps each cell type to its
total cell count. `count_matrix` is the integer matrix used by every
inference function.

The `cluster_filters` keyword drops noisy or off-target clusters before
constructing the count matrix.

### Wide-form (already-pivoted CSV)

A CSV whose columns are cell types and whose rows are clonotypes is
ingested with the single-argument method:

```julia
clono_info, cluster_names, cluster_sizes, count_matrix =
    import_count_matrix("data/Clone_counts_HDM.csv")
```

In this case `clono_info` and `cluster_sizes` are `nothing`.

## Running OU-MCMC inference

Frequencies diffuse along the tree under an Ornstein–Uhlenbeck process —
see [`OUContinuousModel`](@ref) for the full specification.

```julia
plot_init, init_tree, trees, LLs, models, root_ps, upd =
    tree_inference(
        OUContinuousModel(n_samples = 200, burn_in = 2_000, sample_interval = 50),
        cluster_names, count_matrix;
        eqmu = 1.5, eqtheta = 0.1, v = 1.0, d = 0.5, g = 0.5,
    )
```

The return tuple is:

- `plot_init` — `Plots.Plot` of the starting tree (handy for sanity checks),
- `init_tree` — `FelNode` actually used to start MCMC,
- `trees` — `Vector{FelNode}` of posterior topology samples,
- `LLs` — log-likelihoods at the sampled points,
- `models` — sampled `(θ, v, μ)` triples,
- `root_ps` — sampled root state distributions,
- `upd` — the `Update` object holding acceptance ratios and proposal stats.

The keyword arguments to `tree_inference`:

- `eqmu`, `eqtheta`, `v` — initial OU parameters (equilibrium mean,
  mean-reversion strength, variance).
- `d`, `g` — pseudo-counts added to zero entries before the digamma /
  trigamma transforms used to seed leaf Gaussian likelihoods.
- `newt` — an optional pre-built starting tree (otherwise a random tree
  is sampled internally).

## A complete example

The [Worked example](example_workflow.md) page walks through the full
OU-MCMC pipeline (HIPSTR consensus, posterior cloud, UMAP overlay,
animated GIF) on the simulated dataset that ships with the package in
`examples/data/`. The same pipeline also lives as a runnable Jupyter
notebook at [`examples/usage_example.ipynb`](https://github.com/MurrellGroup/Phylotrajectories.jl/blob/main/examples/usage_example.ipynb).
