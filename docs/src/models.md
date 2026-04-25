```@meta
CurrentModule = Phylotrajectories
```

# Model & parameters

`tree_inference` is driven by [`OUContinuousModel`](@ref) — a generative
model in which clonotype log-frequencies diffuse along the tree under an
Ornstein–Uhlenbeck (OU) process: a Brownian motion that mean-reverts
towards an equilibrium. The MCMC sampler updates topology, branch
lengths, the root state, and the OU parameters in turn.

`tree_inference(::OUContinuousModel, cluster_names, cluster_clono_matrix; …)`
expects rows of `cluster_clono_matrix` (a `Matrix{Int64}`) to be cell
types and columns to be clonotypes; `cluster_names::Vector{String}` are
the row labels.

## OUContinuousModel

```@docs
OUContinuousModel
```

Key knobs:

- `tree_warmup_cycles` — number of MCMC iterations with **branch
  lengths and OU parameters frozen** at the start of the run. Useful
  for letting the topology relax first.
- `update::OUContinuousUpdate` — controls the mix of NNI, branch length,
  root, and OU-parameter moves; carries the prior/proposal samplers
  ([`OUEqmuSampler`], [`OUThetaSampler`],
  [`OUGaussianStateSample`](@ref), `BranchlengthSampler`).

`tree_inference(::OUContinuousModel, cluster_names, count_matrix; …)`
additionally accepts:

- `eqmu`, `eqtheta`, `v` — initial OU parameter values (equilibrium
  mean, mean-reversion strength, variance).
- `d`, `g` — pseudo-counts added to zero entries before the digamma /
  trigamma transforms used to seed leaf Gaussian likelihoods.
- `newt` — a pre-built starting tree to bypass the random initialisation.

## Sampler primitives

```@docs
FrequencySampler
RootAcceptanceRatio
GaussianStateSample
OUGaussianStateSample
```
