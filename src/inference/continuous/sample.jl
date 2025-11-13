include("IndependentGaussiansPartition.jl")
include("IndependentBrownianMotion.jl")

"""
    FrequencySampler(proposal::Distribution{Univariate,Continuous})

A type that allows you to specify an additive proposal function. It also holds the acceptance ratio acc_ratio (acc_ratio[1] stores the number of accepts, and acc_ratio[2] stores the number of rejects).
"""
mutable struct FrequencySampler
    acc_ratio::Tuple{Float64, Int64, Int64}
    proposal::Distribution{Univariate,Continuous}
    FrequencySampler(proposal) = new((0.0, 0, 0), proposal)
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
mutable struct RootAcceptanceRatio

Holds tuples of `(ratio, total, #acceptances)` in the fields `position` and `state`, where `ratio::Float64` is the acceptance ratio, `total::Int64` is the total number of proposals, and `#acceptances::Int64` is the number of acceptances.
"""
mutable struct RootAcceptanceRatio
    position::Tuple{Float64, Int64, Int64} #(ratio, total, #acceptances)
    state::Tuple{Float64, Int64, Int64}
    RootAcceptanceRatio() = new((0.0, 0, 0), (0.0, 0, 0))
end

#=
function Base.show(io::IO, r::RootAcceptanceRatio)
    println(io, """\n
Position
    Ratio:   $(r.position[1])
    Total:   $(r.position[2])
    Accepts: $(r.position[3])
State
    Ratio:   $(r.state[1])
    Total:   $(r.state[2])
    Accepts: $(r.state[3])""")
end
=#
Base.show(io::IO, r::RootAcceptanceRatio) = print(io, "position=$(r.position), state=$(r.state)")

"""
# Summary
`struct GaussianStateSample <: MolecularEvolution.UniformRootPositionSample`

Implements the metropolis algorithm for the global Gaussian parameters [μ, ν], where the root state of each clonotype is ~ N(μ, exp(ν)), and updates the root position by `MolecularEvolution.UniformRootPositionSample` if `position = true`. It also holds the acceptance ratio `acc_ratio` in an `RootAcceptanceRatio` struct.
# Constructor
    GaussianStateSample(proposal::ContinuousMultivariateDistribution, prior::ContinuousMultivariateDistribution, radius::Float64, consecutive::Int64; position::Bool = true)

Allows you to specify multivariate proposal and prior distributions for [μ, ν]. `consecutive` is the number of consecutive updates of the root (state *and/or* position) per MCMC iteration. `radius` (∈ [0,1]) scales tree's total branch length to determine local proposal radius for root position.
"""
mutable struct GaussianStateSample{T0, T1<:ContinuousMultivariateDistribution, T2<:ContinuousMultivariateDistribution} <: MolecularEvolution.UniformRootPositionSample
    acc_ratio::RootAcceptanceRatio
    proposal::T1
    prior::T2
    temp_partition::IndependentGaussiansPartition
    radius::Float64
    consecutive::Int64
    parity::Bool #don't mind me, used to track acceptance ratio for position/state
    function GaussianStateSample(
        proposal::T1,
        prior::T2,
        radius::Float64,
        consecutive::Int64;
        position::Bool = true,
    ) where {T1<:ContinuousMultivariateDistribution,T2<:ContinuousMultivariateDistribution}
        @assert length(proposal) == length(prior) == 2 "Proposal and prior must have exactly 2 dimensions"
        @assert 0 <= radius <= 1 "Radius must be in [0, 1]"
        new{position,T1,T2}(RootAcceptanceRatio(), proposal, prior, IndependentGaussiansPartition(0), radius, consecutive, false)
    end
end

const parity_map = Dict(false => :position, true => :state)

Base.length(root_sample::GaussianStateSample) = root_sample.consecutive
MolecularEvolution.radius(root_sample::GaussianStateSample, total_bl::Real) = root_sample.radius * total_bl

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

MolecularEvolution.proposal(modifier::GaussianStateSample{false}, curr_value::@NamedTuple{root::FelNode, dist_above_node::Float64}) = curr_value

function MolecularEvolution.apply_decision(modifier::GaussianStateSample, accept::Bool)
    ratio, total, acc = getproperty(modifier.acc_ratio, parity_map[modifier.parity])
    total += 1
    if accept
        acc += 1
    end
    ratio = acc / total
    setproperty!(modifier.acc_ratio, parity_map[modifier.parity], (ratio, total, acc))
    modifier.parity = !modifier.parity
end

mutable struct MeanDriftSampler{
    T1<:ContinuousUnivariateDistribution,
    T2<:ContinuousUnivariateDistribution,
} <: ModelsUpdate
    acc_ratio::Tuple{Float64, Int64, Int64}
    mean_drift_proposal::T1
    mean_drift_prior::T2
    var_drift::Float64
    function MeanDriftSampler(
        mean_drift_proposal::T1,
        mean_drift_prior::T2,
        var_drift::Float64,
    ) where {T1<:ContinuousUnivariateDistribution, T2<:ContinuousUnivariateDistribution}
        new{T1, T2}((0.0, 0, 0), mean_drift_proposal, mean_drift_prior, var_drift)
    end
end

MolecularEvolution.tr(::MeanDriftSampler, x::IndependentBrownianMotion{Float64}) = x.mean_drifts
MolecularEvolution.invtr(modifier::MeanDriftSampler, x::Float64) =
    IndependentBrownianMotion(x, modifier.var_drift)

MolecularEvolution.proposal(modifier::MeanDriftSampler, curr_value::Float64) =
    curr_value + rand(modifier.mean_drift_proposal)
MolecularEvolution.log_prior(modifier::MeanDriftSampler, x::Float64) =
    logpdf(modifier.mean_drift_prior, x)

# MolecularEvolution.check(modifier::MeanDriftSampler) = true

function (update::MeanDriftSampler)(
    tree::FelNode,
    models::IndependentBrownianMotion{Float64};
    partition_list = 1:length(tree.message),
)
    models = metropolis_step(x->log_likelihood!(tree, x), update, models)

    log_likelihood!(tree, models) #refresh partitions throughout the tree
    return models
end


"""
# Summary
`struct ContinuousUpdate <: MolecularEvolution.AbstractUpdate`

Updates the leaf frequencies, phylogenetic tree, root state and position, and mean drift of the Brownian motion, with metropolis steps.
# Constructor
    ContinuousUpdate(; <keyword arguments>)

# Keyword Arguments
- `branchlength_sampler::MolecularEvolution.BranchlengthSampler=Phylotrajectories.default_branchlength_sampler()`: the proposal and prior distributions for branch length updates in MCMC.
- `frequency_sampler::FrequencySampler=FrequencySampler(Normal())`: the proposal distribution for frequency updates in MCMC.
- `root_sampler::GaussianStateSample=GaussianStateSample(MvNormal(zeros(2), Diagonal([0.1, 0.1])), MvNormal(zeros(2), Diagonal([1.0, 0.1])), 1e-2, 1)`: the proposal and prior distributions for root updates in MCMC.
- `mean_drift_sampler::MeanDriftSampler=MeanDriftSampler(Normal(), Normal(-0.3, 0.5), 1.0)`: the proposal and prior distributions for mean drift updates in MCMC, and the variance drift of the Brownian motion.
- `models::Int=1`: the number of consecutive models updates per MCMC iteration.
- `refresh::Bool=false`: whether to refresh the messages in tree between update operations to ensure message consistency (for debugging purposes).

!!! note
    To disable the sampling of the mean drift, one can set `models=0`.

!!! note
    `GaussianStateSample` also updates the root position by default. See [`GaussianStateSample`](@ref) for more details.
"""
struct ContinuousUpdate <: MolecularEvolution.AbstractUpdate
    bayes_update::MolecularEvolution.StandardUpdate
    temp_messages::Vector{Vector{Partition}}
    refresh::Bool

    function ContinuousUpdate(;
        branchlength_sampler = default_branchlength_sampler(),
        root_sampler = GaussianStateSample(MvNormal(zeros(2), Diagonal([0.1, 0.1])), MvNormal(zeros(2), Diagonal([1.0, 0.1])), 1e-2, 1),
        mean_drift_sampler = MeanDriftSampler(Normal(), Normal(-0.3, 0.5), 1.0),
        models = 1,
        refresh = false,
    )
        new(BayesUpdate(root = 1, models = models, refresh = refresh, branchlength_sampler = branchlength_sampler, root_sampler = root_sampler, models_sampler = mean_drift_sampler), Vector{Vector{Partition}}(), refresh)
    end
end

function (update::ContinuousUpdate)(tree::FelNode, models; partition_list = 1:length(tree.message))
    update.refresh && refresh!(tree, models)
    return update.bayes_update(tree, models, partition_list = partition_list)
end 

MolecularEvolution.collapse_models(::ContinuousUpdate, x::IndependentBrownianMotion) = x.mean_drifts


#### OU Process #######

include("OUpartition.jl")
include("OUprocess.jl")

mutable struct OURootAcceptanceRatio
    position::Tuple{Float64, Int64, Int64}
    state::Tuple{Float64, Int64, Int64}
    OURootAcceptanceRatio() = new((0.0, 0, 0), (0.0, 0, 0))
end
Base.show(io::IO, r::OURootAcceptanceRatio) = print(io, "position=$(r.position), state=$(r.state)")

"""
# Summary
`struct OUGaussianStateSample <: MolecularEvolution.UniformRootPositionSample`

Implements the metropolis algorithm for the global Gaussian parameters [μ, ν], where the root state of each clonotype is ~ N(μ, exp(ν)), and updates the root position by `MolecularEvolution.UniformRootPositionSample` if `position = true`. It also holds the acceptance ratio `acc_ratio` in an `OURootAcceptanceRatio` struct.
# Constructor
    OUGaussianStateSample(proposal::ContinuousMultivariateDistribution, prior::ContinuousMultivariateDistribution, radius::Float64, consecutive::Int64; position::Bool = true)

Allows you to specify multivariate proposal and prior distributions for [μ, ν]. `consecutive` is the number of consecutive updates of the root (state *and/or* position) per MCMC iteration. `radius` (∈ [0,1]) scales tree's total branch length to determine local proposal radius for root position.
"""
mutable struct OUGaussianStateSample{T0, T1<:ContinuousMultivariateDistribution, T2<:ContinuousMultivariateDistribution} <: MolecularEvolution.UniformRootPositionSample
    acc_ratio::OURootAcceptanceRatio
    proposal::T1
    prior::T2
    temp_partition::FBGaussianPartition
    radius::Float64
    consecutive::Int64
    parity::Bool
    function OUGaussianStateSample(
        proposal::T1,
        prior::T2,
        radius::Float64,
        consecutive::Int64;
        position::Bool = true,
    ) where {T1<:ContinuousMultivariateDistribution,T2<:ContinuousMultivariateDistribution}
        @assert length(proposal) == length(prior) == 2 "Proposal and prior must have exactly 2 dimensions"
        @assert 0 <= radius <= 1 "Radius must be in [0, 1]"
        new{position,T1,T2}(OURootAcceptanceRatio(), proposal, prior, FBGaussianPartition(0), radius, consecutive, false)
    end
end

Base.length(root_sample::OUGaussianStateSample) = root_sample.consecutive
MolecularEvolution.radius(root_sample::OUGaussianStateSample, total_bl::Real) = root_sample.radius * total_bl

function set_idg!(dest::FBGaussianPartition, mean::Float64, var::Float64)
    dest.part.mu .= mean
    dest.part.var .= var
    dest.part.log_norm_const .= 0.0
end

gaussian_params(idg::FBGaussianPartition) = idg[1][1:2]

function MolecularEvolution.tr(modifier::OUGaussianStateSample, curr_value::Vector{<:Partition})
    modifier.temp_partition = curr_value[1]
    return tr_gaussian_params(gaussian_params(curr_value[1]))
end

function MolecularEvolution.invtr(modifier::OUGaussianStateSample, tr_curr_value::Vector{Float64})
    set_idg!(modifier.temp_partition, invtr_gaussian_params(tr_curr_value)...)
    return [modifier.temp_partition]
end

MolecularEvolution.proposal(modifier::OUGaussianStateSample, curr_value::Vector{Float64}) = 
    curr_value .+ rand(modifier.proposal)

MolecularEvolution.log_prior(modifier::OUGaussianStateSample, curr_value::Vector{Float64}) = 
    logpdf(modifier.prior, curr_value)

MolecularEvolution.proposal(modifier::OUGaussianStateSample{false}, curr_value::@NamedTuple{root::FelNode, dist_above_node::Float64}) = curr_value

function MolecularEvolution.apply_decision(modifier::OUGaussianStateSample, accept::Bool)
    ratio, total, acc = getproperty(modifier.acc_ratio, parity_map[modifier.parity])
    total += 1
    if accept
        acc += 1
    end
    ratio = acc / total
    setproperty!(modifier.acc_ratio, parity_map[modifier.parity], (ratio, total, acc))
    modifier.parity = !modifier.parity
end


mutable struct OUThetaSampler{
    T1<:ContinuousUnivariateDistribution,
    T2<:ContinuousUnivariateDistribution,
} <: ModelsUpdate
    acc_ratio::Tuple{Float64, Int64, Int64}
    logtheta_proposal::T1
    logtheta_prior::T2
    μ::Float64
    v::Float64
    function OUThetaSampler(
        logtheta_proposal::T1,
        logtheta_prior::T2,
        μ::Float64,
        v::Float64,
    ) where {T1<:ContinuousUnivariateDistribution, T2<:ContinuousUnivariateDistribution}
        new{T1, T2}((0.0, 0, 0), logtheta_proposal, logtheta_prior, μ, v)
    end
end

# Transform: extract log-variance from the model
MolecularEvolution.tr(::OUThetaSampler, x::OrnsteinUhlenbeckModel{Float64}) = log(x.process.θ)
# Inverse transform: construct model from log-variance
MolecularEvolution.invtr(modifier::OUThetaSampler, logtheta::Float64) =
    OrnsteinUhlenbeckModel(ForwardBackward.OrnsteinUhlenbeck(modifier.μ, modifier.v, exp(logtheta)))

# Propose new log-variance
MolecularEvolution.proposal(modifier::OUThetaSampler, curr_logtheta::Float64) =
    curr_logtheta + rand(modifier.logtheta_proposal)
# Log-prior on log-variance
MolecularEvolution.log_prior(modifier::OUThetaSampler, logtheta::Float64) =
    logpdf(modifier.logtheta_prior, logtheta)

# MolecularEvolution.check(modifier::OUThetaSampler) = true

function (update::OUThetaSampler)(
    tree::FelNode,
    models::OrnsteinUhlenbeckModel{Float64};
    partition_list = 1:length(tree.message),
)
    models = metropolis_step(x->log_likelihood!(tree, x), update, models)
    log_likelihood!(tree, models)
    return models
end


mutable struct OUEqmuSampler{
    T1<:ContinuousUnivariateDistribution,
    T2<:ContinuousUnivariateDistribution,
} <: ModelsUpdate
    acc_ratio::Tuple{Float64, Int64, Int64}
    eqmu_proposal::T1
    eqmu_prior::T2
    v::Float64
    θ::Float64
    function OUEqmuSampler(
        eqmu_proposal::T1,
        eqmu_prior::T2,
        v::Float64,
        θ::Float64
    ) where {T1<:ContinuousUnivariateDistribution, T2<:ContinuousUnivariateDistribution}
        new{T1, T2}((0.0, 0, 0), eqmu_proposal, eqmu_prior, v, θ)
    end
end

MolecularEvolution.tr(::OUEqmuSampler, x::OrnsteinUhlenbeckModel{Float64}) = x.process.μ
# Inverse transform: construct model from log-variance
MolecularEvolution.invtr(modifier::OUEqmuSampler, eqmu::Float64) =
    OrnsteinUhlenbeckModel(ForwardBackward.OrnsteinUhlenbeck(eqmu, modifier.v, modifier.θ))

MolecularEvolution.proposal(modifier::OUEqmuSampler, curr_value::Float64) =
    curr_value + rand(modifier.eqmu_proposal)
MolecularEvolution.log_prior(modifier::OUEqmuSampler, eqmu::Float64) =
    logpdf(modifier.eqmu_prior, eqmu)

# MolecularEvolution.check(modifier::OUEqmuSampler) = true

function (update::OUEqmuSampler)(
    tree::FelNode,
    models::OrnsteinUhlenbeckModel{Float64};
    partition_list = 1:length(tree.message),
)
    models = metropolis_step(x->log_likelihood!(tree, x), update, models)

    log_likelihood!(tree, models) #refresh partitions throughout the tree
    return models
end

struct CompositeModelsUpdate <: ModelsUpdate
    eqmu_sampler::ModelsUpdate
    theta_sampler::ModelsUpdate
end

function set_theta!(sampler::OUEqmuSampler, θ::Float64)
    sampler.θ = θ
end

function set_eqmu!(sampler::OUThetaSampler, eqmu::Float64)
    sampler.μ = eqmu
end

function (update::CompositeModelsUpdate)(tree::FelNode, models::OrnsteinUhlenbeckModel;
                                         partition_list=1:length(tree.message))
    # Update theta and variance to current model state
    set_theta!(update.eqmu_sampler, models.process.θ)
    set_eqmu!(update.theta_sampler, models.process.μ)

    # Run samplers
    models = update.eqmu_sampler(tree, models, partition_list=partition_list)
    set_eqmu!(update.theta_sampler, models.process.μ)

    models = update.theta_sampler(tree, models, partition_list=partition_list)
    set_theta!(update.eqmu_sampler, models.process.θ)

    return models
end

struct OUContinuousUpdate <: MolecularEvolution.AbstractUpdate
    bayes_update::MolecularEvolution.StandardUpdate
    temp_messages::Vector{Vector{Partition}}
    refresh::Bool
    function OUContinuousUpdate(;
        branchlength_sampler = default_branchlength_sampler(),
        root_sampler = OUGaussianStateSample(MvNormal(zeros(2), Diagonal([0.1, 0.1])), MvNormal(zeros(2), Diagonal([1.0, 0.1])), 1e-2, 1),
        ou_eqmu_sampler = OUEqmuSampler(Normal(), Normal(1.5, 1.0), 1.0, 0.1),
        ou_theta_sampler = OUThetaSampler(Normal(), Normal(), 1.5, 0.1),
        models = 1,
        refresh = true,
    )
        composite_sampler = CompositeModelsUpdate(ou_eqmu_sampler, ou_theta_sampler)

        new(BayesUpdate(root = 1, 
                        models = models, 
                        refresh = refresh, 
                        branchlength_sampler = branchlength_sampler, 
                        root_sampler = root_sampler, 
                        models_sampler = composite_sampler
                        ),
            Vector{Vector{Partition}}(), refresh)
        end
end

function (update::OUContinuousUpdate)(tree::FelNode, models; partition_list = 1:length(tree.message))
    update.refresh && refresh!(tree, models)
    return update.bayes_update(tree, models, partition_list = partition_list)
end

MolecularEvolution.collapse_models(::OUContinuousUpdate, x::OrnsteinUhlenbeckModel) = [x.process.θ, x.process.v, x.process.μ]