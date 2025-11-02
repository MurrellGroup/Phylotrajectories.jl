########### OrnsteinUhlenbeck Process ############

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