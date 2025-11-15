module Phylotrajectories

using CSV,
    DataFrames, MolecularEvolution, StatsBase, Distributions, Phylo, Plots, LinearAlgebra, ForwardBackward, SpecialFunctions

include("inference/inference.jl")
include("importing.jl")
include("simulations.jl")
include("recombination.jl")

include("viz/plotting.jl")

export OU_MCMC_tree_inference, tree_inference, import_count_matrix, sim_count_matrix, recombine
export IndependentBrownianMotion, IndependentGaussiansPartition, FBGaussianPartition, OrnsteinUhlenbeckModel
export FrequencySampler, RootAcceptanceRatio, GaussianStateSample, MeanDriftSampler, ContinuousUpdate
export OUVarianceSampler, OUContinuousUpdate, OUGaussianStateSample, OUThetaSampler, OUEqmuSampler
export DiscreteModel, ContinuousModel, OUContinuousModel

export PlotTreeOnUmap, PlotTreeOnUmapNoAnimShadow, PlotTreeOnUmapNoAnim

end
