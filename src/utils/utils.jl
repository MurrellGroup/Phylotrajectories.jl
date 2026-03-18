# helper (kept semantically the same as yours)
function collect_marginal_dicts(marg_dict, nsize::Int)
    matrix_states = Matrix{Float64}(undef, nsize, 0)
    # iterate in a deterministic order: by object id of keys
    for key in sort!(collect(keys(marg_dict)); by = x -> objectid(x))
        means = marg_dict[key][1].part.mu
        vars  = marg_dict[key][1].part.var
        states = exp.(means .+ (vars ./ 2))   # log-normal moment
        matrix_states = hcat(matrix_states, states)
    end
    return matrix_states
end


"""
    run_ou_and_build_clone_matrix(
        trees,
        models_tmp,            # iterable of (theta, v, mu) or similar; expects x[1]=θ, x[2]=v, x[3]=μ
        count_matrix::AbstractMatrix,
        cluster_names::AbstractVector{<:AbstractString},
        digamma_value::Real,
        trigamma_value::Real;
        outfile::Union{Nothing,String} = nothing,
        leaf_order::Symbol = :reverse
    ) -> NamedTuple

Builds a HIPSTR tree from `trees`, initializes messages from `count_matrix`,
sets leaf names from `cluster_names`, estimates OU params as medians from `models_tmp`,
computes marginals, and returns the clone matrix (genes × clones), transposed as in your code.

Returns:
`(hip, ou_model, marg_state_dict, clone_matrix, start_ll)`

If `outfile` is a string, writes the clone matrix as CSV.
"""

function run_ou_and_build_clone_matrix(
    trees,
    models_tmp,
    count_matrix::AbstractMatrix,
    cluster_names::AbstractVector{<:AbstractString},
    digamma_value::Real,
    trigamma_value::Real,
    outfile::String
)
    # --- basic checks
    @assert !isempty(trees) "trees is empty."
    @assert size(count_matrix, 1) == length(cluster_names) "Rows of count_matrix must match length of cluster_names."

    # 1) Tidy trees and build HIPSTR container
    ladderize!.(trees)
    hip = HIPSTR(trees)

    # 2) OU params as medians from models_tmp (expects x[1]=θ, x[2]=v, x[3]=μ as in your snippet)
    mu    = median(getindex.(models_tmp, 3))
    theta = median(getindex.(models_tmp, 1))
    v     = median(getindex.(models_tmp, 2))

    # 3) Init message template on HIPSTR object
    message_template = [FBGaussianPartition(size(count_matrix, 2))]
    internal_message_init!(hip, message_template)

    # 4) Digamma/Trigamma “pseudo-count” handling for zeros
    cluster_clono_matrix_digamma  = ifelse.(count_matrix .== 0, count_matrix .+ digamma_value,  count_matrix)
    cluster_clono_matrix_trigamma = ifelse.(count_matrix .== 0, count_matrix .+ trigamma_value, count_matrix)

    # 5) Set leaf names and update leaf partitions from counts
    leaves = getleaflist(hip)
    @assert length(leaves) == length(cluster_names) "Number of leaves in HIPSTR must match cluster_names."

    for (i, n) in enumerate(reverse(leaves))
        n.name = cluster_names[i]

        digamma_values  = digamma.(cluster_clono_matrix_digamma[i, :])
        trigamma_values = trigamma.(cluster_clono_matrix_trigamma[i, :])

        update_leaf_partition!(n,
            digamma_values,
            trigamma_values,
            Float64.(count_matrix[i, :]))
    end

    # 6) Aggregate parent messages across input trees via medians (means/vars/consts)
    #    (kept identical to your approach for reproducibility)
    median_means  = median.(eachrow(hcat(getfield.(getfield.(first.(getfield.(trees, :parent_message)), :part), :mu)...)))
    median_vars   = median.(eachrow(hcat(getfield.(getfield.(first.(getfield.(trees, :parent_message)), :part), :var)...)))
    median_consts = median.(eachrow(hcat(getfield.(getfield.(first.(getfield.(trees, :parent_message)), :part), :log_norm_const)...)))

    # Set the "parent message" on the current HIPSTR tree
    update_leaf_partition!(hip, median_means, median_vars, median_consts)

    # 7) Build OU model and get starting log-likelihood
    ou_model = Phylotrajectories.OrnsteinUhlenbeckModel(
        ForwardBackward.OrnsteinUhlenbeck(mu, v, theta)
    )
    
    start_ll = log_likelihood!(hip, ou_model)
    # println("Starting LL: ", start_ll)

    # 8) Marginals → clone matrix
    marg_state_dict = marginal_state_dict(hip, ou_model)

    clone_matrix = collect_marginal_dicts(marg_state_dict, size(count_matrix, 2))'
    # shape: (#marginals) × (#clones); same as your code (note the final transpose)

    CSV.write(outfile, DataFrame(clone_matrix, :auto))
end


tuple_sort(t::Tuple{UInt64,UInt64}) = ifelse(t[1] < t[2], t, (t[2], t[1]))

function node_hash_split(node, hash_container, node_container, name2hash; push_leaves = false)
    if isleafnode(node)
        if push_leaves
            push!(hash_container, name2hash[node.name])
            push!(node_container, node)
        end
        return name2hash[node.name]
    else
        child_hashes = [
            node_hash_split(nc, hash_container, node_container, name2hash, push_leaves = push_leaves) for
            nc in node.children
        ]
        #merge child hashes with xor as we go up the tree
        first_hash = child_hashes[1]
        for ch in child_hashes[2:end]
            first_hash = xor(first_hash, ch)
        end
        push!(hash_container, first_hash)
        push!(node_container, node)
        return first_hash
    end
end

function get_node_hashes_rooted(newt; push_leaves = false)
    leafnames = [n.name for n in getleaflist(newt)]
    leafhashes = hash.(leafnames)
    all_names_hash = xor(leafhashes...)
    name2hash = Dict(zip(leafnames, leafhashes))
    hash_container = UInt64[]
    node_container = FelNode[]
    #This puts things in the containers
    node_hash_split(newt, hash_container, node_container, name2hash, push_leaves = push_leaves)
    #This makes a hash that matches everything except the given node
    other_hash = xor.(hash_container, all_names_hash)
    #Sort these, to make comparisons order invariant, which makes the comparison rooting invariant
    #Only sort internal node hashes, and don't use other_hash for leaves (to avoid collisions)
    isleafposition = isleafnode.(node_container)
    sensitive_tuple_sort(t::Tuple{Bool, Tuple{UInt64, UInt64}}) = ifelse(t[1], (t[2][1], t[2][1]), tuple_sort(t[2]))
    #Consider making this sort an option, and then we can have a rooted comparison and an unrooted one
    sorted_hash_pairs = collect(zip(isleafposition, zip(hash_container, other_hash)))
    return sorted_hash_pairs, node_container
end

#returns nodes in the query that don't have matching splits in the reference
function tree_diff_rooted(query, reference)
    newt_hc, newt_nc = get_node_hashes_rooted(query)
    n_hc, n_nc = get_node_hashes_rooted(reference)
    hashset = Set(n_hc)
    changed_nodes = newt_nc[[!(n in hashset) for n in newt_hc]]
    return changed_nodes
end

function getname(node::FelNode)
    return node.name
end

function VectorOfDistances(tree)

    distmat, node_dic = tree2distances(tree)
    leaflist = getleaflist(tree)
    sort!(leaflist, by = x -> x.name)
    order = [node_dic[leaf] for leaf in leaflist]
    distmat = distmat[order, order]
    node_distances = distmat[triu(trues(size(distmat)), 1)] 

    root_distances, node_dic = root2tip_distances(tree)
    order = [node_dic[leaf] for leaf in leaflist]
    root_distances = root_distances[order]

    return vcat(node_distances, root_distances)
end

### _______________ Metrics ____________________
function SimComparison(newtree, oldtree)
    x = VectorOfDistances(newtree);
    y = VectorOfDistances(oldtree);

    pearson_cor = Statistics.cor(Float64.(x),Float64.(y))
    spearman_cor = corspearman(Float64.(x),Float64.(y))

    node_differences = tree_diff_rooted(oldtree, newtree)
    append!(node_differences, tree_diff_rooted(newtree, oldtree))

    return pearson_cor, spearman_cor, length(node_differences)
end
