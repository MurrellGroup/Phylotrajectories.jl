module Phylotrajectories

using CSV,
    DataFrames, MolecularEvolution, StatsBase, Distributions, Phylo, Plots, LinearAlgebra

include("inference/discrete/discrete.jl")
include("inference/continuous/continuous.jl")
include("importing.jl")
include("simulations.jl")
include("recombination.jl")

export tree_inference, import_count_matrix, sim_count_matrix, recombine
export IndependentBrownianMotion, IndependentGaussiansPartition
end
