const CANONICAL_JUMP = 0.3

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
    tree_inference(cluster_names::Vector{String}, cluster_clono_matrix::Matrix{Int64}; <keyword arguments>)

Find the Maximum Likelihood tree for a count matrix `cluster_clono_matrix`.
Returns the ML tree, CTMC model, the discretized states over frequencies, and the LL of the tree

# Arguments
...
"""
function tree_inference(
    cluster_names::Vector{String},
    cluster_clono_matrix::Matrix{Int64};
    jump = CANONICAL_JUMP,
    a = 1.0,
    b = 1.0,
    Ne = 1.0,
    rate = 10.0,
    start_branch_length = 0.1,
    max_cycles = 10,
    n_random_trees = 1,
)
    @assert size(cluster_names, 1) == size(cluster_clono_matrix, 1)
    lowest_average = minimum(mean(cluster_clono_matrix, dims = 2)) # mean per leaves -> minimum
    states = exp.([log(lowest_average)-0.5:jump:log(maximum(cluster_clono_matrix))+1;])
    println("Number of states ", length(states))

    model = DiagonalizedCTMC(Q(length(states), a, b) .* CANONICAL_JUMP / jump)
    message_template =
        [CustomDiscretePartition(length(states), size(cluster_clono_matrix)[2])]

    LLs = Vector{Float64}()
    ML_newt = FelNode()
    MLL = -Inf
    for _ = 1:n_random_trees
        #Random starting tree
        newt = sim_tree(size(cluster_clono_matrix)[1], Ne, rate)
        internal_message_init!(newt, message_template)

        #Set the leaf names from the imported count matrix, and init the partitions based on the counts there
        for (i, n) in enumerate(reverse(getleaflist(newt)))
            n.name = cluster_names[i]
            poisson_partition!(n.message[1], states, cluster_clono_matrix[i, :])
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


        #Set the parent message.
        #This sets the Q states that correspond to low counts to have some probability mass at the root.
        #This is like an inductive bias that the root will have "naive" unexpanded cells.
        #We might want to do something a bit more elegant, and actually try and learn these quentities, but the signal might not be there.
        thresh_ind = findfirst(states .> 1)
        newt.parent_message[1].state[1:thresh_ind, :] .= 1 / sum(thresh_ind)

        println("Starting LL: ", log_likelihood!(newt, model))

        #Optimize the tree topology:
        @time for i = 1:max_cycles #Needs a stopping condition check. I should add a "topology only" search to MolecularEvolution.jl...
            felsenstein_down!(newt, model)
            nni_optim!(newt, model)
            # println("LL: ", log_likelihood!(newt, model))
        end

        #Polish topology and branch lengths:
        @time tree_polish!(newt, model, verbose = 0, tol = 10^-6)
        ladderize!(newt)
        LL = log_likelihood!(newt, model)
        if LL > MLL
            MLL = LL
            ML_newt = newt
        end
        push!(LLs, LL)
    end

    return ML_newt, model, states, MLL, LLs
end
