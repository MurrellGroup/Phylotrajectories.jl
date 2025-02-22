module Phylotrajectories

using CSV,
    DataFrames, MolecularEvolution, StatsBase, Distributions, Phylo, Plots, LinearAlgebra

include("inference/inference.jl")
include("importing.jl")
include("simulations.jl")
include("recombination.jl")

export tree_inference, import_count_matrix, sim_count_matrix, recombine
export IndependentBrownianMotion, IndependentGaussiansPartition
export FrequencySampler, GaussianStateSample, MeanDriftSampler, ContinuousUpdate
export DiscreteModel, ContinuousModel
end
