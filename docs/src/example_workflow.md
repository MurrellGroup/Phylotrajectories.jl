```@meta
CurrentModule = Phylotrajectories
```

# Worked example: OU-MCMC + UMAP overlay on simulated data

This page walks through the analysis pipeline that produced the
*DogAllergen* figures in the paper, but on a small simulated dataset that
ships with the package. The same workflow is also available as a
runnable Jupyter notebook at
[`examples/usage_example.ipynb`](https://github.com/MurrellGroup/Phylotrajectories.jl/blob/main/examples/usage_example.ipynb).

The simulated data (`examples/data/simulated_clone_data.tsv` and
`examples/data/simulated_UMAP_coords.tsv`) are produced by
[`examples/generate_simulated_data.jl`](https://github.com/MurrellGroup/Phylotrajectories.jl/blob/main/examples/generate_simulated_data.jl)
using [`sim_count_matrix`](@ref).

## Set-up

```julia
using Phylotrajectories, MolecularEvolution, ForwardBackward
using CSV, DataFrames, Distributions, Statistics, LinearAlgebra
using Plots, Plots.PlotMeasures, ColorSchemes, Phylo, Random

Random.seed!(1234)
data_dir = joinpath(pkgdir(Phylotrajectories), "examples", "data")
path     = mktempdir()              # where outputs go
```

## 1. Load the simulated counts

```julia
clono_info, cluster_names, cluster_sizes, count_matrix = import_count_matrix(
    joinpath(data_dir, "simulated_clone_data.tsv"),
    :Clonotype, :cell_types, :TRB_cdr3aa,
)
```

## 2. Define the OU MCMC model

```julia
function make_ou_model(; burn_in = 2_000, sample_interval = 50,
                          n_samples = 200, tree_warmup_cycles = 100)
    OUContinuousModel(
        update = OUContinuousUpdate(
            nni = 1, branchlength = 1, root = 1, models = 1,
            branchlength_sampler = BranchlengthSampler(Normal(0, 0.1), Normal(-1, 1)),
            root_sampler = OUGaussianStateSample(
                MvNormal(zeros(2), Diagonal([0.01, 0.01])),
                MvNormal(zeros(2), Diagonal([1.0, 0.1])),
                1e-1, 1),
            ou_eqmu_sampler  = OUEqmuSampler(Normal(0.0, 2.0), Normal(1.5, 1.0), 1.0, 0.1),
            ou_theta_sampler = OUThetaSampler(Normal(0, 0.5), Normal(-1, 1), 1.5, 1.0),
        ),
        start_branch_length = 0.1,
        tree_warmup_cycles  = tree_warmup_cycles,
        burn_in             = burn_in,
        sample_interval     = sample_interval,
        n_samples           = n_samples,
    )
end

plot_init, init_tree, trees, LLs, models, root_ps, upd =
    tree_inference(make_ou_model(), cluster_names, count_matrix;
                   eqmu = 1.5, eqtheta = 0.1, v = 1.0,
                   d = 0.5, g = 0.5)
```

## 3. HIPSTR consensus + posterior cloud

```julia
ladderize!.(trees)
hip, node2logcred, node2support = HIPSTR(trees; getcred = true, getsupport = true)
mltpl = plot_multiple_trees(trees, hip; line_width = 0.075)
```

## 4. Build the clone × node frequency matrix

```julia
run_ou_and_build_clone_matrix(
    trees, models, count_matrix, cluster_names,
    0.5, 0.5,
    joinpath(path, "OU_sim_clone_matrix.csv"),
)
```

## 5. Project the inferred tree onto a UMAP

The package ships with three plotting helpers for overlaying a tree on a
2-D embedding:

| Helper                                  | What it draws                                        |
|-----------------------------------------|------------------------------------------------------|
| [`PlotTreeOnUmapNoAnimShadow`](@ref)    | Faint branches (used for an MCMC posterior cloud).   |
| [`PlotTreeOnUmapNoAnim`](@ref)          | Solid branches + nodes (used for the consensus tree). |
| [`PlotTreeOnUmap`](@ref)                | Animated walk; returns a `Plots.@animate` object.    |

The intended pattern is to:

1. Build a base scatter of the cells in UMAP space.
2. For every posterior tree, copy it, pin the leaves to their cell-type
   centroids, run `felsenstein!` under a 2-D Brownian motion, and add the
   tree as a shadow (`PlotTreeOnUmapNoAnimShadow`).
3. Do the same once more for the HIPSTR consensus, drawing it with
   solid lines (`PlotTreeOnUmapNoAnim`).

The full code is reproduced in the notebook —
[`examples/usage_example.ipynb`](https://github.com/MurrellGroup/Phylotrajectories.jl/blob/main/examples/usage_example.ipynb).
