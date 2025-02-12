include("IndependentGaussiansPartition.jl")

"""
    FrequencySampler(proposal::Distribution{Univariate,Continuous})

A type that allows you to specify an additive proposal function. It also holds the acceptance ratio acc_ratio (acc_ratio[1] stores the number of accepts, and acc_ratio[2] stores the number of rejects).
"""
struct FrequencySampler
    acc_ratio::Array{Int64,1}
    proposal::Distribution{Univariate,Continuous}
    FrequencySampler(proposal) = new(zeros(Int64, 2), proposal)
end

MolecularEvolution.proposal(modifier::FrequencySampler, curr_value::Array{Float64,1}) =
    curr_value + rand(modifier.proposal, length(curr_value))
MolecularEvolution.log_prior(modifier::FrequencySampler, x) = 0.0 #Constant, since, prior comes from root
function MolecularEvolution.apply_decision(modifier::FrequencySampler, accepts::BitArray)
    for b in accepts
        MolecularEvolution.apply_decision(modifier, b)
    end
end
#This should idiomatically be achieved by a loop or broadcasting metropolis_step over multiple functions and Float64s
#but the driving factor for why I specialize metropolis_step is
#it's convenient for LL to return Array{Float64,1}
function MolecularEvolution.metropolis_step(
    LL::Function,
    modifier::FrequencySampler,
    curr_value::Array{Float64,1},
)
    prop = MolecularEvolution.proposal(modifier, curr_value)
    accept_proposal =
        rand(length(curr_value)) .<=
        exp.(
            LL(prop) .+ MolecularEvolution.log_prior(modifier, prop) .- LL(curr_value) .-
            MolecularEvolution.log_prior(modifier, curr_value),
        )
    MolecularEvolution.apply_decision(modifier, accept_proposal)
    return ifelse.(accept_proposal, prop, curr_value)
end

"""
# Summary
`struct GaussianStateSample <: MolecularEvolution.UniformRootPositionSample`

Implements the metropolis algorithm for the global Gaussian parameters [μ, ν], where the root state of each clonotype is ~ N(μ, exp(ν)), and updates the root position by `MolecularEvolution.UniformRootPositionSample`.
# Constructor
    GaussianStateSample(proposal::ContinuousMultivariateDistribution, prior::ContinuousMultivariateDistribution, temp_partition::IndependentGaussiansPartition)

Allows you to specify multivariate proposal and prior distributions for [μ, ν]. It also holds the acceptance ratio `acc_ratio` (`acc_ratio[1]` stores the number of accepts, and `acc_ratio[2]` stores the number of rejects).
"""
mutable struct GaussianStateSample{T1, T2} <: MolecularEvolution.UniformRootPositionSample where {T1,T2 <: ContinuousMultivariateDistribution}
    acc_ratio::Array{Int64,1}
    proposal::T1
    prior::T2
    temp_partition::IndependentGaussiansPartition
    function GaussianStateSample(
        proposal::T1,
        prior::T2,
    ) where {T1<:ContinuousMultivariateDistribution,T2<:ContinuousMultivariateDistribution}
        @assert length(proposal) == length(prior) == 2 "Proposal and prior must have exactly 2 dimensions"
        new{T1,T2}(zeros(Int64, 2), proposal, prior, IndependentGaussiansPartition(0))
    end
end

function set_idg!(dest::IndependentGaussiansPartition, mean::Float64, var::Float64)
    dest.means .= mean
    dest.vars .= var
    dest.norm_consts .= 0.0
end

gaussian_params(idg::IndependentGaussiansPartition) = idg[1][1:2]
# Transform mean and variance to mean and log-variance 
function tr_gaussian_params(curr_values)
    curr_values .|> [identity, log]
end

# Inverse transform mean and log-variance to mean and variance
function invtr_gaussian_params(tr_curr_values)
    tr_curr_values .|> [identity, exp]
end

function MolecularEvolution.tr(modifier::GaussianStateSample, curr_value::Vector{<:Partition})
    modifier.temp_partition = curr_value[1]
    return tr_gaussian_params(gaussian_params(curr_value[1]))
end

function MolecularEvolution.invtr(modifier::GaussianStateSample, tr_curr_value::Vector{Float64})
    set_idg!(modifier.temp_partition, invtr_gaussian_params(tr_curr_value)...)
    return [modifier.temp_partition]
end

MolecularEvolution.proposal(modifier::GaussianStateSample, curr_value::Vector{Float64}) = 
    curr_value .+ rand(modifier.proposal)

MolecularEvolution.log_prior(modifier::GaussianStateSample, curr_value::Vector{Float64}) = 
    logpdf(modifier.prior, curr_value)

"""
# Summary
`struct ContinuousUpdate <: MolecularEvolution.AbstractUpdate`

Updates the leaf frequencies, phylogenetic tree, and root distribution with metropolis steps.
# Constructor
    ContinuousUpdate(; <keyword arguments>)

# Keyword Arguments
- `branchlength_sampler::MolecularEvolution.BranchlengthSampler=DEFAULT_BRANCHLENGTH_SAMPLER`: the proposal and prior distributions for branch length updates in MCMC.
- `frequency_sampler::FrequencySampler=FrequencySampler(Normal())`: the proposal distribution for frequency updates in MCMC.
- `root_sampler::GaussianStateSample=GaussianStateSample(MvNormal(zeros(2), Diagonal([0.1, 0.1])), MvNormal(zeros(2), Diagonal([1.0, 0.1])))`: the proposal and prior distributions for root updates in MCMC.

!!! note
    `GaussianStateSample` also updates the root position. See [`GaussianStateSample`](@ref) for more details.
"""
struct ContinuousUpdate <: MolecularEvolution.AbstractUpdate
    bayes_update::MolecularEvolution.StandardUpdate
    frequency_sampler::FrequencySampler
    temp_messages::Vector{Vector{Partition}}

    function ContinuousUpdate(;
        branchlength_sampler = DEFAULT_BRANCHLENGTH_SAMPLER,
        frequency_sampler = FrequencySampler(Normal()),
        root_sampler = GaussianStateSample(MvNormal(zeros(2), Diagonal([0.1, 0.1])), MvNormal(zeros(2), Diagonal([1.0, 0.1]))),
    )
        new(BayesUpdate(root = 1, branchlength_sampler = branchlength_sampler, root_sampler = root_sampler), frequency_sampler, Vector{Vector{Partition}}())
    end
end

function (update::ContinuousUpdate)(tree::FelNode, models; partition_list = 1:length(tree.message))
    sample_leafs!(update.temp_messages, tree, x -> models, update.frequency_sampler)
    return update.bayes_update(tree, models, partition_list = partition_list)
end