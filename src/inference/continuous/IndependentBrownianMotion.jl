#Model behavior
mutable struct IndependentBrownianMotion <: MolecularEvolution.ContinuousStateModel
    mean_drifts::Array{Float64,1}
    var_drifts::Array{Float64,1}
    function IndependentBrownianMotion(
        mean_drifts::Array{Float64,1},
        var_drifts::Array{Float64,1},
    )
        @assert length(mean_drifts) == length(var_drifts)
        new(mean_drifts, var_drifts)
    end
    # Use constant mean and var drift for all independent gaussians
    # Works because for example [1, 2] .* [1] = [1, 2]
    # Defined internally to avoid @assert
    function IndependentBrownianMotion(mean_drift::Float64, var_drift::Float64)
        new([mean_drift], [var_drift])
    end
end

function MolecularEvolution.backward!(
    dest::IndependentGaussiansPartition,
    source::IndependentGaussiansPartition,
    model::IndependentBrownianMotion,
    node::FelNode,
)
    dest.means .= source.means .- node.branchlength .* model.mean_drifts
    dest.vars .= source.vars .+ node.branchlength .* model.var_drifts
    dest.norm_consts .= source.norm_consts
end

function MolecularEvolution.forward!(
    dest::IndependentGaussiansPartition,
    source::IndependentGaussiansPartition,
    model::IndependentBrownianMotion,
    node::FelNode,
)
    dest.means .= source.means .+ node.branchlength .* model.mean_drifts
    dest.vars .= source.vars .+ node.branchlength .* model.var_drifts
    dest.norm_consts .= source.norm_consts
end

#If you want to use a root prior, you can set these values explicitly.
function MolecularEvolution.eq_freq_from_template(
    model::IndependentBrownianMotion,
    partition_template::IndependentGaussiansPartition,
)
    out_partition = copy_partition(partition_template)
    out_partition.means .= 0.0
    out_partition.vars .= Inf
    out_partition.norm_consts .= 0.0
    return out_partition
end
