# Phylotrajectories

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://MurrellGroup.github.io/Phylotrajectories.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://MurrellGroup.github.io/Phylotrajectories.jl/dev/)
[![Build Status](https://github.com/MurrellGroup/Phylotrajectories.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/MurrellGroup/Phylotrajectories.jl/actions/workflows/CI.yml?query=branch%3Amain)


## Importing data

To import count matrices from your data files, you can use the `import_count_matrix` function. Below are examples of how to import data from different file formats.

TODO: show where these files can be found.

### Importing from a TSV file with filters

```julia
using Phylotrajectories

clono_info, cluster_names, count_matrix = import_count_matrix(
    "data/clone_data_HDM.tsv",
    :Clonotype,
    :cell_types,
    :TRB_cdr3aa,
    cluster_filters = ["Proliferating"],
)
```
In this file, `:Clonotype`, `:cell_types` and `:TRB_cdr3aa` are what the clonotype column, cluster column, and cdr3 column are called respectively. We can use the `cluster_filters` keyword to ignore data from specific clusters (in this case, the `"Proliferating"` cluster).

### Importing from a CSV file
```julia
clono_info, cluster_names, count_matrix = import_count_matrix("data/Clone_counts_HDM.csv")
```
> [!NOTE]  
> This time, `clono_info` is `nothing`.

## Performing inference

You can perform inference using either a discrete or continuous model:

> [!NOTE]  
> `DiscreteModel` supports both inferring the Maximum Likelihood tree and approximating a posterior tree space via MCMC, whereas `ContinuousModel` only supports the latter.
```julia
# Using a discrete model
model = DiscreteModel()
newtree, model, states, LL, LLs = tree_inference(model, cluster_names, count_matrix)

# Using a continuous model
model = ContinuousModel()
newtree, model, trees, LLs = tree_inference(model, cluster_names, count_matrix)
```

See [`DiscreteModel`](@ref) and [`ContinuousModel`](@ref) for available model parameters.

TODO: Add tree-plot fig.