module Phylotrajectories

using CSV, DataFrames, MolecularEvolution, StatsBase, Distributions, Phylo, Plots

include("inference.jl")
include("importing.jl")
#TODO: Simulations?

export tree_inference, import_count_matrix
end
