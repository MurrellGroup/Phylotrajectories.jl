```@meta
CurrentModule = Phylotrajectories
```

# Worked example: OU-MCMC on a simulated count matrix

This page walks through the OU-MCMC pipeline on the small simulated
dataset that ships with the package
(`examples/data/simulated_clone_data.csv`).  The same workflow is
available as a runnable notebook at
[`examples/usage_example.ipynb`](https://github.com/MurrellGroup/Phylotrajectories.jl/blob/main/examples/usage_example.ipynb).

## Set-up

```julia
using Phylotrajectories, MolecularEvolution, ForwardBackward
using CSV, DataFrames, Distributions, Statistics, StatsBase, LinearAlgebra
using Plots, Plots.PlotMeasures, ColorSchemes, Phylo, Random

Random.seed!(1234)
data_dir = joinpath(pkgdir(Phylotrajectories), "examples", "data")
path     = mktempdir()
```

## 1. Read the simulated clone count matrix

The bundled `simulated_clone_data.csv` is a wide-form CSV: rows are
clonotypes, columns are cell-type subsets (`Subset1`…`Subset6`).
The single-argument form of [`import_count_matrix`](@ref) returns a
4-tuple `(clono_info, cluster_names, cluster_sizes, count_matrix)` —
with `clono_info` and `cluster_sizes` set to `nothing` for wide-form
inputs, since they don't carry per-cell metadata.

```julia
_, cluster_names, _, count_matrix = import_count_matrix(
    joinpath(data_dir, "simulated_clone_data.csv"),
)
```

## 2. Subsample the clonotypes

The full matrix has ~1 900 clonotypes — more than this example needs.
We retain 750 of them so the OU MCMC finishes in minutes:

```julia
sampled_indices = sample(1:size(count_matrix, 2), 750; replace = false)
count_matrix = count_matrix[:, sampled_indices]

@show cluster_names                       # cell-type labels (matrix rows)
@show sum(count_matrix; dims = 2) |> vec  # per-subset totals — should all be > 0
```

A row that sums to zero would make the OU likelihood degenerate, so this
is also the right place to reject-sample if needed.

## 3. Define the OU MCMC model

```julia
digamma_value  = 0.5
trigamma_value = 0.5

eqmu    = 1.5
eqtheta = 0.1
v       = 1.0

function make_ou_model(; burn_in = 5_000, sample_interval = 100,
                          n_samples = 100, tree_warmup_cycles = 100)
    OUContinuousModel(
        update = OUContinuousUpdate(
            nni = 1, branchlength = 1, root = 1, models = 1,
            branchlength_sampler = BranchlengthSampler(Normal(0, 0.1), Normal(-1, 1)),
            root_sampler = OUGaussianStateSample(
                MvNormal(zeros(2), Diagonal([0.01, 0.01])),
                MvNormal(zeros(2), Diagonal([1.0, 0.1])),
                1e-1, 1),
            ou_eqmu_sampler  = OUEqmuSampler(Normal(0.0, 2.0), Normal(1.5, 1.0), 1.0, 0.1),
            ou_theta_sampler = OUThetaSampler(Normal(0, 0.1), Normal(-1, 1), 1.5, 1.0),
        ),
        start_branch_length = 0.1,
        tree_warmup_cycles  = tree_warmup_cycles,
        burn_in             = burn_in,
        sample_interval     = sample_interval,
        n_samples           = n_samples,
    )
end
```

`digamma_value` and `trigamma_value` are pseudo-counts added to zero
entries before the digamma / trigamma transforms used to seed the leaf
Gaussian likelihoods. `eqmu`, `eqtheta`, `v` are the initial values of
the OU process (equilibrium mean, mean-reversion strength, variance).
`tree_warmup_cycles` runs the topology-only update at the start so the
tree relaxes before the parameter updates kick in.

## 4. Run inference

```julia
plot_init, init_tree, trees, LLs, models, root_ps, upd =
    tree_inference(make_ou_model(),
                   cluster_names, count_matrix;
                   eqmu = eqmu, eqtheta = eqtheta, v = v,
                   d = digamma_value, g = trigamma_value)
```

## 5. HIPSTR consensus + posterior cloud

`HIPSTR` builds a maximum-credibility tree from the posterior samples
and returns the per-node support values for the credibility plot.

```julia
ladderize!.(trees)
hip, node2logcred, node2support = HIPSTR(trees; getcred = true, getsupport = true)
mltpl = plot_multiple_trees(trees, hip; line_width = 1.0)
```

## 6. Diagnostics & credibility plot

The notebook saves a 6-panel diagnostic dashboard (LL trace, OU `θ`,
variance, `EqMu` acceptance ratios, root-state cloud, MCMC posterior
cloud overlaid on HIPSTR) and a HIPSTR tree with annotated per-node
posterior support.  The full code is in
[`examples/usage_example.ipynb`](https://github.com/MurrellGroup/Phylotrajectories.jl/blob/main/examples/usage_example.ipynb).
