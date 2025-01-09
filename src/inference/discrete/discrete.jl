function Q(numstates, a, b)
    Qmat = zeros(numstates, numstates)
    Qmat[1, 1] = -b
    Qmat[1, 2] = b
    Qmat[end, end] = -a
    Qmat[end, end-1] = a
    for i = 2:size(Qmat)[1]-1
        Qmat[i, i] = -(a + b)
        Qmat[i, i+1] = b
        Qmat[i, i-1] = a
    end
    return Qmat
end

function poisson_partition!(
    dest::DiscretePartition,
    rate_list::Array{Float64,1},
    observed_counts::Array{Int64,1},
)
    for i = 1:length(observed_counts)
        dest.state[:, i] = [pdf(Poisson(s), observed_counts[i]) for s in rate_list]
    end
end

"""
    tree_inference(model::DiscreteModel, cluster_names::Vector{String}, cluster_clono_matrix::Matrix{Int64})

Find the Maximum Likelihood tree for a count matrix `cluster_clono_matrix` using a discrete state model.
Returns the ML tree, CTMC model, discretized states of frequencies, optimized trees, and the final LLs of the initial trees.

See [`DiscreteModel`](@ref) for model specification and parameters.
"""
function tree_inference(
    model::DiscreteModel,
    cluster_names::Vector{String},
    cluster_clono_matrix::Matrix{Int64},
)
    @assert size(cluster_names, 1) == size(cluster_clono_matrix, 1)
    lowest_average = minimum(mean(cluster_clono_matrix, dims = 2))
    states =
        exp.([log(lowest_average)-0.5:model.jump:log(maximum(cluster_clono_matrix))+1;])
    println("Number of states ", length(states))

    ctmc_model = DiagonalizedCTMC(
        Q(length(states), model.a, model.b) .* (CANONICAL_JUMP / (model.jump^2)),
    )
    message_template =
        [CustomDiscretePartition(length(states), size(cluster_clono_matrix)[2])]
    
    if model.ML
        #--------------------------------
        #ML
        #--------------------------------
        LLs = Vector{Float64}()
        ML_newt = FelNode()
        MLL = -Inf
        trees = Vector{FelNode}()
        for _ = 1:model.n_random_trees
            newt = sim_tree(size(cluster_clono_matrix)[1], model.Ne, model.sample_rate)
            push!(trees, newt)
            internal_message_init!(newt, message_template)

            for (i, n) in enumerate(reverse(getleaflist(newt)))
                n.name = cluster_names[i]
                poisson_partition!(n.message[1], states, cluster_clono_matrix[i, :])
            end

            i = 1
            for n in getnodelist(newt)
                n.nodeindex = i
                if !MolecularEvolution.isroot(n)
                    n.branchlength = model.start_branch_length
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

            thresh_ind = findfirst(states .> 1)
            newt.parent_message[1].state[1:thresh_ind, :] .= 1 / sum(thresh_ind)

            println("Starting LL: ", log_likelihood!(newt, ctmc_model))

            @time for i = 1:model.max_cycles
                felsenstein_down!(newt, ctmc_model)
                nni_optim!(newt, ctmc_model)
            end

            @time tree_polish!(newt, ctmc_model, verbose = 0, tol = 10^-6)
            ladderize!(newt)
            LL = log_likelihood!(newt, ctmc_model)
            if LL > MLL
                MLL = LL
                ML_newt = newt
            end
            push!(LLs, LL)
        end
        return ML_newt, ctmc_model, states, trees, LLs
    else
        #--------------------------------
        #MCMC
        #--------------------------------
        newt = sim_tree(size(cluster_clono_matrix)[1], model.Ne, model.sample_rate)
        internal_message_init!(newt, message_template)

        for (i, n) in enumerate(reverse(getleaflist(newt)))
            n.name = cluster_names[i]
            poisson_partition!(n.message[1], states, cluster_clono_matrix[i, :])
        end

        i = 1
        for n in getnodelist(newt)
            n.nodeindex = i
            if !MolecularEvolution.isroot(n)
                n.branchlength = model.start_branch_length
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

        thresh_ind = findfirst(states .> 1)
        newt.parent_message[1].state[1:thresh_ind, :] .= 1 / sum(thresh_ind)

        println("Starting LL: ", log_likelihood!(newt, ctmc_model))
        trees, LLs = metropolis_sample(
            newt,
            [ctmc_model],
            model.n_samples,
            bl_sampler = model.branchlength_sampler,
            burn_in = model.burn_in,
            sample_interval = model.sample_interval,
            collect_LLs = true,
        )
        return newt, ctmc_model, states, trees, LLs
    end
end
