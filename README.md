# Phylotrajectories

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://nossleinad.github.io/Phylotrajectories.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://nossleinad.github.io/Phylotrajectories.jl/dev/)
[![Build Status](https://github.com/nossleinad/Phylotrajectories.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/nossleinad/Phylotrajectories.jl/actions/workflows/CI.yml?query=branch%3Amain)


## Importing Data

To import count matrices from your data files, you can use the `import_count_matrix` function. Below are examples of how to import data from different file formats.

### Importing from a TSV file with filters

```julia
using Phylotrajectories, MolecularEvolution

clono_info, cluster_names, count_matrix = import_count_matrix(
    "data/clone_data_HDM.tsv",
    :Clonotype,
    :cell_types,
    :TRB_cdr3aa,
    cluster_filters = ["Proliferating"],
)
```

### Importing from a CSV file
```julia
clono_info, cluster_names, count_matrix = import_count_matrix("data/Clone_counts_HDM.csv")
```
> [!NOTE]  
> This time, clono_info is nothing.

## Performing inference
```julia
newtree, model, states, LL = tree_inference(cluster_names, count_matrix)
```

## Simulating a count matrix
Coming soon...