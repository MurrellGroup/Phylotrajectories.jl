include("IndependentGaussiansPartition.jl")
include("IndependentBrownianMotion.jl")

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

Implements the metropolis algorithm for the global Gaussian parameters [μ, ν], where the root state of each clonotype is ~ N(μ, exp(ν)), and updates the root position by `MolecularEvolution.UniformRootPositionSample`. It also holds the acceptance ratio `acc_ratio` (`acc_ratio[1]` stores the number of accepts, and `acc_ratio[2]` stores the number of rejects).
# Constructor
    GaussianStateSample(proposal::ContinuousMultivariateDistribution, prior::ContinuousMultivariateDistribution, consecutive::Int64)

Allows you to specify multivariate proposal and prior distributions for [μ, ν]. `consecutive` is the number of consecutive updates of the root (state *and* position) per MCMC iteration.
"""
mutable struct GaussianStateSample{T1, T2} <: MolecularEvolution.UniformRootPositionSample where {T1,T2 <: ContinuousMultivariateDistribution}
    acc_ratio::Array{Int64,1}
    proposal::T1
    prior::T2
    temp_partition::IndependentGaussiansPartition
    consecutive::Int64
    function GaussianStateSample(
        proposal::T1,
        prior::T2,
        consecutive::Int64,
    ) where {T1<:ContinuousMultivariateDistribution,T2<:ContinuousMultivariateDistribution}
        @assert length(proposal) == length(prior) == 2 "Proposal and prior must have exactly 2 dimensions"
        new{T1,T2}(zeros(Int64, 2), proposal, prior, IndependentGaussiansPartition(0), consecutive)
    end
end

Base.length(root_sample::GaussianStateSample) = root_sample.consecutive

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

struct MeanDriftSampler{
    T1<:ContinuousUnivariateDistribution,
    T2<:ContinuousUnivariateDistribution,
} <: ModelsUpdate
    acc_ratio::Vector{Int}
    mean_drift_proposal::T1
    mean_drift_prior::T2
    var_drift::Float64
    function MeanDriftSampler(
        mean_drift_proposal::T1,
        mean_drift_prior::T2,
        var_drift::Float64,
    ) where {T1<:ContinuousUnivariateDistribution, T2<:ContinuousUnivariateDistribution}
        new{T1, T2}([0, 0], mean_drift_proposal, mean_drift_prior, var_drift)
    end
end

MolecularEvolution.tr(::MeanDriftSampler, x::IndependentBrownianMotion{Float64}) = x.mean_drifts
MolecularEvolution.invtr(modifier::MeanDriftSampler, x::Float64) =
    IndependentBrownianMotion(x, modifier.var_drift)

MolecularEvolution.proposal(modifier::MeanDriftSampler, curr_value::Float64) =
    curr_value + rand(modifier.mean_drift_proposal)
MolecularEvolution.log_prior(modifier::MeanDriftSampler, x::Float64) =
    logpdf(modifier.mean_drift_prior, x)

function (update::MeanDriftSampler)(
    tree::FelNode,
    models::IndependentBrownianMotion{Float64};
    partition_list = 1:length(tree.message),
)
    metropolis_step(update, models) do x::IndependentBrownianMotion{Float64}
        log_likelihood!(tree, x)
    end
end


"""
# Summary
`struct ContinuousUpdate <: MolecularEvolution.AbstractUpdate`

Updates the leaf frequencies, phylogenetic tree, root state and position, and mean drift of the Brownian motion, with metropolis steps.
# Constructor
    ContinuousUpdate(; <keyword arguments>)

# Keyword Arguments
- `branchlength_sampler::MolecularEvolution.BranchlengthSampler=DEFAULT_BRANCHLENGTH_SAMPLER`: the proposal and prior distributions for branch length updates in MCMC.
- `frequency_sampler::FrequencySampler=FrequencySampler(Normal())`: the proposal distribution for frequency updates in MCMC.
- `root_sampler::GaussianStateSample=GaussianStateSample(MvNormal(zeros(2), Diagonal([0.1, 0.1])), MvNormal(zeros(2), Diagonal([1.0, 0.1])), 1)`: the proposal and prior distributions for root updates in MCMC.
- `mean_drift_sampler::MeanDriftSampler=MeanDriftSampler(Normal(), Normal(-0.3, 0.5), 1.0)`: the proposal and prior distributions for mean drift updates in MCMC, and the variance drift of the Brownian motion.
- `models::Int=1`: the number of consecutive models updates per MCMC iteration.

!!! note
    To disable the sampling of the mean drift, one can set `models=0`.

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
        root_sampler = GaussianStateSample(MvNormal(zeros(2), Diagonal([0.1, 0.1])), MvNormal(zeros(2), Diagonal([1.0, 0.1])), 1),
        mean_drift_sampler = MeanDriftSampler(Normal(0.0, 0.2), Normal(-0.3, 0.5), 1.0),
        models = 1
    )
        new(BayesUpdate(root = 1, models = models, branchlength_sampler = branchlength_sampler, root_sampler = root_sampler, models_sampler = mean_drift_sampler), frequency_sampler, Vector{Vector{Partition}}())
    end
end

function (update::ContinuousUpdate)(tree::FelNode, models; partition_list = 1:length(tree.message))
    sample_leafs!(update.temp_messages, tree, x -> [models], update.frequency_sampler)
    return update.bayes_update(tree, models, partition_list = partition_list)
end 

MolecularEvolution.collapse_models(::ContinuousUpdate, x::IndependentBrownianMotion) = x.mean_drifts