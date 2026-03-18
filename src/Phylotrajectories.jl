module Phylotrajectories

using CSV,
    DataFrames, MolecularEvolution, StatsBase, Distributions, Phylo, Plots, LinearAlgebra, ForwardBackward, SpecialFunctions, Statistics

include("inference/inference.jl")
include("importing.jl")
include("simulations.jl")
include("recombination.jl")

include("viz/plotting.jl")
include("utils/utils.jl")

export tree_inference, import_count_matrix, sim_count_matrix, recombine
export IndependentBrownianMotion, IndependentGaussiansPartition, FBGaussianPartition, OrnsteinUhlenbeckModel, LogExponential
export FrequencySampler, RootAcceptanceRatio, GaussianStateSample, MeanDriftSampler, ContinuousUpdate
export OUVarianceSampler, OUContinuousUpdate, OUGaussianStateSample, OUThetaSampler, OUEqmuSampler
export DiscreteModel, ContinuousModel, OUContinuousModel

export PlotTreeOnUmap, PlotTreeOnUmapNoAnimShadow, PlotTreeOnUmapNoAnim
export run_ou_and_build_clone_matrix, SimComparison

end
