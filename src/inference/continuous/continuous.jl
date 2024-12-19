include("IndependentGaussiansPartition.jl")
include("IndependentBrownianMotion.jl")
include("sample.jl")

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

function continuous_tree_inference(
    cluster_names::Vector{String},
    cluster_clono_matrix::Matrix{Int64};
    mean_drift = 1.0,
    Ne = 1.0,
    sample_rate = 10.0,
    start_branch_length = 0.1,
    samples = 10,
    consecutive_root_samples = 10,
)
    @assert size(cluster_names, 1) == size(cluster_clono_matrix, 1)
    model = IndependentBrownianMotion(mean_drift, 1.0)
    message_template = [IndependentGaussiansPartition(size(cluster_clono_matrix)[2])]


    #Random starting tree
    newt = sim_tree(size(cluster_clono_matrix)[1], Ne, sample_rate)
    internal_message_init!(newt, message_template)

    #Set the leaf names from the imported count matrix, and init the partitions based on the counts there
    for (i, n) in enumerate(reverse(getleaflist(newt)))
        n.name = cluster_names[i]
        #obs2partition!(n.message[1], log.(cluster_clono_matrix[i, :] .+ 1)) #Start with MLE (avoid log(0))
        n.message[1].vars .= 0.0
        poisson_partition!(n.message[1], cluster_clono_matrix[i, :])
    end

    i = 1
    for n in getnodelist(newt)
        n.nodeindex = i
        if !MolecularEvolution.isroot(n)
            n.branchlength = start_branch_length
        end

        i += 1
    end

    ladderize!(newt)
    phylo_tree = get_phylo_tree(newt)
    plot(
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
    )

    #We let the prior over root frequencies be Normal(0, 1)

    println("Starting LL: ", log_likelihood!(newt, model))

    temp_messages = [copy_message(newt.message)]
    #Run the MCMC with custom leaf sampling
    trees, LLs, = metropolis_sample(
        newt,
        [model],
        samples,
        burn_in = 2000,
        collect_LLs = true,
    ) do tree, models
        sample_leafs!(temp_messages, tree, x -> models)
        for i = 1:consecutive_root_samples
            sample_root_distribution!(temp_messages[1], tree)
        end
        nni_optim!(tree, x -> models, selection_rule = softmax_sampler)
        branchlength_optim!(tree, x -> models)
    end


    return newt, model, trees, LLs
end

#Have temp_messages as input so we can reuse between chain iterations
function sample_leafs!(
    temp_messages::Vector{Vector{T}},
    tree::FelNode,
    models;
    partition_list = 1:length(tree.message),
    sampler::FrequencySampler = FrequencySampler(Normal()),
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

#Assumes that we've done an up pass
function sample_root_distribution!(
    temp_message::Vector{<:Partition},
    tree::FelNode;
    sampler::GaussianSampler = GaussianSampler(
        MvNormal(zeros(2), Diagonal(ones(2))),
        MvNormal(zeros(2), Diagonal(ones(2))),
    ),
)
    LL(x) = root_LL(x, temp_message, tree)
    gaussian_params = metropolis_step(LL, sampler, collect(tree.parent_message[1][1][1:2]))
    set_idg!(tree.parent_message[1], gaussian_params...)
end

function root_LL(
    gaussian_params::Array{Float64,1},
    temp_message::Vector{<:Partition},
    tree::FelNode,
)
    part = 1
    set_idg!(temp_message[part], gaussian_params...)
    combine!(temp_message[part], tree.message[part])

    return sum(MolecularEvolution.total_LL.(temp_message))
end

function set_idg!(dest::IndependentGaussiansPartition, mean::Float64, var::Float64)
    dest.means .= mean
    dest.vars .= var
    dest.norm_consts .= 0.0
end
