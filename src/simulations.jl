"""
    sim_count_matrix(
        tree::FelNode,
        n_clonotypes::Int,
        n_cells::Int,
        initial_partition::Partition,
        models;
        outpath = "",
    )

Simulate a count matrix for a given tree. Returns the cluster names and a `cell_types`-by-`clonotypes` count matrix, where `cell_types` is the number of leaves on the tree and `clonotypes <= n_clonotypes` is the number of clonotypes whose sampled total count in the tree is non-zero. The count matrix sums to `n_cells`.
`models` can either be a single `BranchModel` (if the model is global for all clonotypes) or an array of `BranchModel`s, if the models are clonotype-wise, or 
a function that takes a clonotype index, and returns a `FelNode -> Vector{<:BranchModel}` function if you need clonotype-wise models to vary from one branch to another.
"""
function sim_count_matrix(
    tree::FelNode,
    n_clonotypes::Int,
    n_cells::Int,
    initial_partition::Partition,
    models;
    outpath = "",
)
    @assert partition2obs(initial_partition) isa Float64
    n_cell_types = length(getleaflist(tree))
    freqs = Matrix{Float64}(undef, n_cell_types, n_clonotypes)
    internal_message_init!(tree, initial_partition)
    for i = 1:n_clonotypes
        model = models(i)
        sample_down!(tree, model)
        freqs[:, i] = [partition2obs(n.message[1]) for n in getleaflist(tree)]
    end
    #Softmax the logits to get probabilities
    expfreqs = exp.(freqs)
    expfreqs ./= sum(expfreqs)

    #Sample from the probabilities to get counts
    flat_draw = rand(Multinomial(n_cells, reshape(expfreqs, :)))

    count_matrix = reshape(flat_draw, n_cell_types, :)

    leaf_names = [n.name for n in getleaflist(tree)]
    no_zeros = count_matrix[:, sum(count_matrix, dims = 1)[:].>0]
    df = DataFrame()
    for i = 1:n_cell_types
        df[!, leaf_names[i]] = no_zeros[i, :]
    end

    if outpath != ""
        CSV.write(outpath * "_count_matrix.csv", df)
        open(outpath * ".tre", "a") do io
            write(io, newick(tree))
        end
    end
    return leaf_names, no_zeros
end

function sim_count_matrix(
    tree::FelNode,
    n_clonotypes::Int,
    n_cells::Int,
    initial_partition::Partition,
    models::Vector{<:BranchModel};
    outpath = "",
)
    @assert length(models) == n_clonotypes
    sim_count_matrix(
        tree,
        n_clonotypes,
        n_cells,
        initial_partition,
        i -> models[i],
        outpath = outpath,
    )
end

sim_count_matrix(
    tree::FelNode,
    n_clonotypes::Int,
    n_cells::Int,
    initial_partition::Partition,
    model::BranchModel;
    outpath = "",
) = sim_count_matrix(
    tree,
    n_clonotypes,
    n_cells,
    initial_partition,
    i -> model;
    outpath = outpath,
)

#This should maybe be in MolecularEvolution.jl
function getnode2rootpath(node::T) where {T<:AbstractTreeNode}
    nodelist = T[node]
    currnode = node
    while !isnothing(currnode.parent)
        currnode = currnode.parent
        push!(nodelist, currnode)
    end
    return nodelist
end
