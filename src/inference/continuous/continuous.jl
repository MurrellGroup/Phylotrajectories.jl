#################### OU MCMC tree inference ####################

function MolecularEvolution.metropolis_sample(
    update!::AbstractUpdate,
    initial_tree::FelNode,
    models,#::Vector{<:BranchModel},
    num_of_samples;
    partition_list = 1:length(initial_tree.message),
    burn_in = 1000,
    sample_interval = 10,
    collect_LLs = false,
    midpoint_rooting = false,
    ladderize = false,
    collect_models = false,
    )

    # The prior over the (log) of the branchlengths should be specified in bl_sampler. 
    # Furthermore, a non-informative/uniform prior is assumed over the tree topolgies (excluding the branchlengths).

    sample_LLs = Float64[]
    samples = FelNode[]
    sample_models = []
    tree = initial_tree
    iterations = burn_in + num_of_samples * sample_interval

    #matrix to collect marginal_state_dict from nodes of sampled trees
    root_params = []

    for i = 1:iterations
        # Updates the tree topolgy and branchlengths.
        tree, models = update!(tree, models, partition_list = partition_list)
        if isnothing(tree)
            break
        end

        if (i - burn_in) % sample_interval == 0 && i > burn_in

            if collect_models
                push!(sample_models, collapse_models(update!, models))
            end

            if collect_LLs
                push!(sample_LLs, log_likelihood!(tree, models, partition_list = partition_list))
            end

            push!(samples, deepcopy(tree))
            push!(root_params, tree.parent_message[1][1][1:2])

        end

    end

    if midpoint_rooting
        for (i, sample) in enumerate(samples)
            node, len = midpoint(sample)
            samples[i] = reroot!(node, dist_above_child = len)
        end
    end

    if ladderize
        for sample in samples
            ladderize!(sample)
        end
    end

    if collect_LLs && collect_models
        return samples, sample_LLs, sample_models, root_params
    elseif collect_LLs && !collect_models
        return samples, sample_LLs
    elseif !collect_LLs && collect_models
        return samples, sample_models
    end


    return samples
end


# Helper function to update leaf nodes
function update_leaf_partition!(
    node::FelNode,
    digamma_values::Vector{Float64},
    trigamma_values::Vector{Float64},
    log_norm_const::Vector{Float64}
)
    # Create new GaussianLikelihood with the updated values
    new_part = ForwardBackward.GaussianLikelihood(
        copy(digamma_values),
        copy(sqrt.(trigamma_values)),
        copy(log_norm_const)
    )
    
    # Update the message with a new FBGaussianPartition
    node.message[1] = FBGaussianPartition(new_part)
end

"""
    tree_inference(model::ContinuousModel, cluster_names::Vector{String}, cluster_clono_matrix::Matrix{Int64})

Returns the initial tree, sampled trees, LLs of the sampled trees, and the mean drift of the sampled models.

See [`ContinuousModel`](@ref) for model specification and parameters.
"""

function OU_MCMC_tree_inference(
    model::OUContinuousModel,
    cluster_names::Vector{String},
    cluster_clono_matrix::Matrix{Int64};
    newt = missing,
    d=0.5, 
    g=0.5,
    eqtheta = 1.0,
    eqmu = 1.5,
    v = 1.0
)
    @assert size(cluster_names, 1) == size(cluster_clono_matrix, 1)
    ou_model = Phylotrajectories.OrnsteinUhlenbeckModel(ForwardBackward.OrnsteinUhlenbeck(eqmu, v, eqtheta))
    message_template = [FBGaussianPartition(size(cluster_clono_matrix)[2])]

    # We add pseudocounts to zero counts to avoid numerical issues with the Poisson likelihood
    if ismissing(newt)
        println("Using random tree.")
        newt = sim_tree(size(cluster_clono_matrix)[1], model.Ne, model.sample_rate)
        cluster_clono_matrix_digamma = ifelse.(cluster_clono_matrix .== 0, cluster_clono_matrix .+ d, cluster_clono_matrix)
        cluster_clono_matrix_trigamma = ifelse.(cluster_clono_matrix .== 0, cluster_clono_matrix .+ g, cluster_clono_matrix)

        internal_message_init!(newt, message_template)

        #Set the leaf names from the imported count matrix, and init the partitions based on the counts there
        for (i, n) in enumerate(reverse(getleaflist(newt)))
            n.name = cluster_names[i]

            digamma_values = digamma.(cluster_clono_matrix_digamma[i, :])
            trigamma_values = trigamma.(cluster_clono_matrix_trigamma[i, :])

            update_leaf_partition!(n, digamma_values, trigamma_values, Float64.(cluster_clono_matrix[i, :]))
        end

        i = 1
        for n in getnodelist(newt)
            n.nodeindex = i
            if !MolecularEvolution.isroot(n)
                n.branchlength =  model.start_branch_length
        end

        i += 1
        end

        #Set the parent message.
        update_leaf_partition!(newt, 
        repeat([0.0]  , size(cluster_clono_matrix)[2]), 
        repeat([1.0], size(cluster_clono_matrix)[2]), 
        repeat([0.], size(cluster_clono_matrix)[2]))
    else
        println("Using existing tree.")
    end

    ladderize!(newt)
    phylo_tree = get_phylo_tree(newt)
    plot_init = plot(
        phylo_tree,
        showtips = true,
        tipfont = 6,
        markersize = 4.0,
        markerstrokewidth = 0,
        margins = 1Plots.cm,
        linewidth = 1.5,
        markercolor = :black,
        size = (500, 500),
        title = "Starting Tree",
    );

    println("Starting LL: ", log_likelihood!(newt, ou_model))

    # push!(model.update.temp_messages, copy_message(newt.message))
    print("Inference")
    trees, LLs, models, root_ps = metropolis_sample(
        model.update,
        newt,
        ou_model,
        model.n_samples,
        burn_in = model.burn_in,
        sample_interval = model.sample_interval,
        collect_LLs = true,
        collect_models = true,
    )

    return plot_init, newt, trees, LLs, models, root_ps, model.update.bayes_update
end


#################### MCMC tree inference #################### 

function poisson_partition!(
    dest::IndependentGaussiansPartition,
    observed_counts::Array{Int64,1},
)
    dest.counts = observed_counts
    poisson_partition!(dest)
end

function poisson_partition!(dest::IndependentGaussiansPartition)
    dest.norm_consts .= logpdf.(Poisson.(exp.(dest.means)), dest.counts)
end

"""
    tree_inference(model::ContinuousModel, cluster_names::Vector{String}, cluster_clono_matrix::Matrix{Int64})

Returns the initial tree, sampled trees, LLs of the sampled trees, and the mean drift of the sampled models.

See [`ContinuousModel`](@ref) for model specification and parameters.
"""
function tree_inference(
    model::ContinuousModel,
    cluster_names::Vector{String},
    cluster_clono_matrix::Matrix{Int64};
    newt = missing
)
    @assert size(cluster_names, 1) == size(cluster_clono_matrix, 1)
    bm_model = IndependentBrownianMotion(model.mean_drift, 1.0)
    message_template = [IndependentGaussiansPartition(size(cluster_clono_matrix)[2])]

    if ismissing(newt)
        println("random tree")
        newt = sim_tree(size(cluster_clono_matrix)[1], model.Ne, model.sample_rate)
        internal_message_init!(newt, message_template)

        for (i, n) in enumerate(reverse(getleaflist(newt)))
            n.name = cluster_names[i]
            n.message[1].vars .= 0.0
            poisson_partition!(n.message[1], cluster_clono_matrix[i, :])
        end

        i = 1
        for n in getnodelist(newt)
            n.nodeindex = i
            if !MolecularEvolution.isroot(n)
                n.branchlength = model.start_branch_length
            end
            i += 1
        end
    else
        println("Using existing tree.")
    end


    ladderize!(newt)
    phylo_tree = get_phylo_tree(newt)
    plot_init = plot(
        phylo_tree,
        showtips = true,
        tipfont = 6,
        markersize = 4.0,
        markerstrokewidth = 0,
        margins = 1Plots.cm,
        linewidth = 1.5,
        markercolor = :black,
        size = (500, 500),
        title = "Starting Tree",
    );

    println("Starting LL: ", log_likelihood!(newt, bm_model))

    push!(model.update.temp_messages, copy_message(newt.message))

    trees, LLs, models, root_ps = metropolis_sample(
        model.update,
        newt,
        bm_model,
        model.n_samples,
        burn_in = model.burn_in,
        sample_interval = model.sample_interval,
        collect_LLs = true,
        collect_models = true,
    )

     return plot_init, newt, trees, LLs, models, root_ps, model.update.bayes_update #this is more of a MolEv concern, but I think we're interested in the Log posterior instead of LL?
end

#The message-preserving algorithm is copy-pasted from branchlength/nni_optim! with some minor tweaks that does the leaf sampling.
#TODO: (perhaps) make an interface for message-preserving traversal where you'd just pass in a function that does the leaf sampling.
function sample_leafs!(
    temp_messages::Vector{Vector{T}},
    tree::FelNode,
    models,
    sampler::FrequencySampler,
    partition_list = 1:length(tree.message),
    traversal = Iterators.reverse,
) where {T<:Partition}
    stack = [(pop!(temp_messages), tree, 1, 1, true, true)]
    while !isempty(stack)
        temp_message, node, ind, lastind, first, down = pop!(stack)
        #We start out with a regular downward pass...
        #(except for some extra bookkeeping to track if node is visited for the first time)
        #-------------------
        if !isleafnode(node)
            if down
                if first
                    model_list = models(node)
                    for part in partition_list
                        forward!(
                            temp_message[part],
                            node.parent_message[part],
                            model_list[part],
                            node,
                        )
                    end
                    #Temp must be constant between iterations for a node during down...
                    child_iter = traversal(1:length(node.children))
                    lastind = Base.first(child_iter) #(which is why we track the last child to be visited during down)
                    push!(stack, (Vector{T}(), node, ind, lastind, false, false)) #... but not up
                    for i in child_iter #Iterative reverse <=> Recursive non-reverse, also optimal for lazysort!??
                        push!(stack, (temp_message, node, i, lastind, false, true))
                    end
                end
                if !first
                    sib_inds = sibling_inds(node.children[ind])
                    for part in partition_list
                        combine!(
                            (node.children[ind]).parent_message[part],
                            [mess[part] for mess in node.child_messages[sib_inds]],
                            true,
                        )
                        combine!(
                            (node.children[ind]).parent_message[part],
                            [temp_message[part]],
                            false,
                        )
                    end
                    #But calling sample_leafs! recursively... (the iterative equivalent)
                    push!(
                        stack,
                        (
                            MolecularEvolution.safepop!(temp_messages, temp_message),
                            node.children[ind],
                            ind,
                            lastind,
                            true,
                            true,
                        ),
                    ) #first + down combination => safepop!
                    ind == lastind && push!(temp_messages, temp_message) #We no longer need constant temp
                end
            end
            if !down
                #Then combine node.child_messages into node.message...
                for part in partition_list
                    combine!(
                        node.message[part],
                        [mess[part] for mess in node.child_messages],
                        true,
                    )
                end
                #But now we need to prop back up to set your parents children message correctly.
                #-------------------
                if !MolecularEvolution.isroot(node)
                    temp_message = pop!(temp_messages)
                    model_list = models(node)
                    #Then we need to set the "message_to_set", which is node.parent.child_messages[but_the_right_one]
                    for part in partition_list
                        backward!(
                            node.parent.child_messages[ind][part],
                            node.message[part],
                            model_list[part],
                            node,
                        )
                    end
                    push!(temp_messages, temp_message)
                end
            end
        else
            #But now we need to optimize the current node, and then prop back up to set your parents children message correctly.
            #-------------------
            model_list = models(node)
            fun = x -> leaf_log_posterior(x, temp_message, node, model_list)
            frequencies = metropolis_step(fun, sampler, node.message[1].means)
            node.message[1].means .= frequencies
            poisson_partition!(node.message[1])
            #Consider checking for improvement, and bailing if none.
            #Then we need to set the "message_to_set", which is node.parent.child_messages[but_the_right_one]
            for part in partition_list
                backward!(
                    node.parent.child_messages[ind][part],
                    node.message[part],
                    model_list[part],
                    node,
                )
            end
            push!(temp_messages, temp_message)
        end
    end
end

function leaf_log_posterior(
    frequencies::Array{Float64,1},
    temp_message::Vector{<:Partition},
    node::FelNode,
    model_list::Vector{<:BranchModel},
)
    part = 1
    node.message[1].means = frequencies #not . since we don't want to overwrite
    poisson_partition!(node.message[1])

    backward!(temp_message[part], node.message[part], model_list[part], node)
    combine!(temp_message[part], node.parent_message[part])

    return copy(MolecularEvolution.site_LLs(temp_message[part]))
end
