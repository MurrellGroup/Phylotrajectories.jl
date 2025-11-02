

############ ForwardBackward GaussianLikelihood ############


#TODO: consider slipping in more @assert for dimension matching

#-----------------Gaussian prop---------------
#Partition behavior
mutable struct FBGaussianPartition <: MolecularEvolution.ContinuousPartition
    part::ForwardBackward.GaussianLikelihood

    # Inner constructor to ensure we always have a valid GaussianLikelihood
    function FBGaussianPartition(part::ForwardBackward.GaussianLikelihood{Float64})
        new(part)
    end
end

# Constructors
function FBGaussianPartition(mu::Vector{Float64}, var::Vector{Float64}, log_norm_const::Vector{Float64})
    @assert length(mu) == length(var) == length(log_norm_const)
    FBGaussianPartition(ForwardBackward.GaussianLikelihood(mu, var, log_norm_const))
end

FBGaussianPartition(mu::Vector{Float64}, var::Vector{Float64}) = 
    FBGaussianPartition(mu, var, zeros(length(mu)))

FBGaussianPartition(n::Int) =
    FBGaussianPartition(zeros(n), ones(n), zeros(n))

Base.length(g::FBGaussianPartition) = length(g.part.mu)
Base.length(part::ForwardBackward.GaussianLikelihood) = length(part.mu)

function MolecularEvolution.states(g::FBGaussianPartition) 
    return 0
end

function Base.getindex(g::FBGaussianPartition, i::Int)
    1 <= i <= length(g) || throw(BoundsError(g, i))
    return g.part.mu[i], g.part.var[i], g.part.log_norm_const[i]
end

function Base.setindex!(
    g::FBGaussianPartition,
    v::Tuple{Float64,Float64,Float64},
    i::Int,
)
    1 <= i <= length(g) || throw(BoundsError(g, i))
    g.part.mu[i] = v[1]
    g.part.var[i] = v[2]
    g.part.log_norm_const[i] = v[3]
end

#Overloading the copy_partition to avoid deepcopy.
function MolecularEvolution.copy_partition(src::FBGaussianPartition)
        return FBGaussianPartition(
        copy(src.part.mu),
        copy(src.part.var),
        copy(src.part.log_norm_const),
    )
end

function MolecularEvolution.copy_partition_to!(
    dest::FBGaussianPartition,
    src::FBGaussianPartition,
)
    dest.part.mu .= src.part.mu
    dest.part.var .= src.part.var
    dest.part.log_norm_const .= src.part.log_norm_const
end

function MolecularEvolution.identity!(dest::FBGaussianPartition)
    dest.part.var .= Inf
    dest.part.log_norm_const .= 0.0
end

function MolecularEvolution.combine!(
    dest::FBGaussianPartition,
    src::FBGaussianPartition,
)
    dest.part = dest.part ⊙ src.part
end

function MolecularEvolution.site_LLs(part::FBGaussianPartition)
    return part.part.log_norm_const
end

function Distributions.logpdf(
    g::FBGaussianPartition,
    x::Array{Float64,1},
)
    result = logpdf.(Normal.(g.part.mu, sqrt.(g.part.var)), x)
    for i = 1:length(g)
        if g.part.var[i] == 0
            error("logpdf not defined for point mass")
        end
    end
    return result
end

#And sampling
function MolecularEvolution.sample_partition!(partition::FBGaussianPartition)
    partition.part.mu .= randn(length(partition)) .* sqrt.(partition.part.var) .+ partition.part.mu
    partition.part.var .= 0.0
    partition.part.log_norm_const .= 0.0
end

#And max
function MolecularEvolution.max_partition!(partition::FBGaussianPartition)
    partition.part.var .= 0.0
    partition.part.log_norm_const .= 0.0
end

function MolecularEvolution.obs2partition!(
    partition::FBGaussianPartition,
    obs::Array{Float64,1},
)
    partition.part.mu .= obs
    partition.part.var .= 0.0
    partition.part.log_norm_const .= 0.0
end

function MolecularEvolution.partition2obs(partition::FBGaussianPartition)
    return partition.part.mu
end
