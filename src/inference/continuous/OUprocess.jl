########### OrnsteinUhlenbeck Process ############

struct LogExponential{T<:Real} <: ContinuousUnivariateDistribution
    λ::T
    function LogExponential(λ::T) where {T<:Real}
        λ > zero(T) || throw(ArgumentError("rate λ must be positive"))
        new{T}(λ)
    end
end
Distributions.params(d::LogExponential) = (d.λ,)
Base.minimum(::LogExponential) = -Inf
Base.maximum(::LogExponential) = Inf
function Distributions.insupport(::LogExponential, z::Real)
    return isfinite(z)
end
function Distributions.logpdf(d::LogExponential, z::Real)
    λ = d.λ
    return log(λ) + z - λ * exp(z)
end
function Distributions.cdf(d::LogExponential, z::Real)
    λ = d.λ
    return 1 - exp(-λ * exp(z))
end
function Distributions.quantile(d::LogExponential, p::Real)
    λ = d.λ
    return log(-log(1 - p) / λ)
end

#Model behavior
mutable struct OrnsteinUhlenbeckModel{T} <: MolecularEvolution.ContinuousStateModel where {T<:Union{Float64, Array{Float64,1}}}
    process::ForwardBackward.OrnsteinUhlenbeck{T}
end



function MolecularEvolution.backward!(
    dest::FBGaussianPartition,
    source::FBGaussianPartition,
    model::OrnsteinUhlenbeckModel,
    node::FelNode,
)
    ForwardBackward.backward!(dest.part, source.part, model.process, node.branchlength)
end

function MolecularEvolution.forward!(
    dest::FBGaussianPartition,
    source::FBGaussianPartition,
    model::OrnsteinUhlenbeckModel,
    node::FelNode,
)
    ForwardBackward.forward!(dest.part, source.part, model.process, node.branchlength)
end