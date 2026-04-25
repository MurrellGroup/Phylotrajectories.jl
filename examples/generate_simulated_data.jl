"""
generate_simulated_data.jl

Generate the simulated dataset that backs `examples/usage_example.ipynb`:

    examples/data/simulated_clone_data.tsv   # one row per cell (Barcode, cell_types, Clonotype, ...)
    examples/data/simulated_UMAP_coords.tsv  # per-cell UMAP coordinates

The data are produced with `Phylotrajectories.sim_count_matrix`, which
diffuses log-frequencies down a random tree under Brownian motion and then
draws multinomial counts.  A "bias" Brownian motion is layered onto a
chosen clonotype/cell-type pair so that the simulated set has visible
clone–phenotype structure (matching `docs/src/simulations.md`).

Run from the package root:

    julia --project=. examples/generate_simulated_data.jl
"""

using Phylotrajectories
using MolecularEvolution
using Random, Distributions, Statistics
using DataFrames, CSV

const SEED = 20260425
Random.seed!(SEED)

# ---------------------------------------------------------------------------
# 1. Simulation knobs
# ---------------------------------------------------------------------------
const CELL_TYPES = [
    "Naive", "Tfh", "Th2", "Th17", "Tregs",
    "Act_circ", "Th_CTL_like", "Act2",
]
const N_CELL_TYPES = length(CELL_TYPES)
const N_CLONOTYPES = 40
const N_CELLS      = 4_000

# ---------------------------------------------------------------------------
# 2. Build a random coalescent-like tree with `N_CELL_TYPES` leaves
# ---------------------------------------------------------------------------
n(t) = (10 * N_CELL_TYPES) / (1 + exp(t - 10))
tree = sim_tree(N_CELL_TYPES, n, N_CELL_TYPES / 5; mutation_rate = 0.05)
ladderize!(tree)

# Replace the random leaf labels with our cell-type names
for (i, leaf) in enumerate(getleaflist(tree))
    leaf.name = CELL_TYPES[i]
end

# ---------------------------------------------------------------------------
# 3. Simulate counts via Brownian motion on the tree
# ---------------------------------------------------------------------------
initial_partition = GaussianPartition(0.0, 0.0)
default_BM_model  = BrownianMotion(-0.3, 1.5)

# Pick a (clonotype, cell-type) pair that should display strong bias
bias_clonotype  = 7
towards_cluster = 3
pos_bias_BM_model = BrownianMotion( 1.5, 1.5)
neg_bias_BM_model = BrownianMotion(-4.0, 1.5)

function bias_model(i)
    if i != bias_clonotype
        return default_BM_model
    end
    bias_branches = Set(Phylotrajectories.getnode2rootpath(getleaflist(tree)[towards_cluster]))
    d = Dict{FelNode,BrownianMotion}()
    for nd in getnodelist(tree)
        d[nd] = nd ∈ bias_branches ? pos_bias_BM_model : neg_bias_BM_model
    end
    return n::FelNode -> [d[n]]
end

cluster_names, count_matrix = sim_count_matrix(
    tree, N_CLONOTYPES, N_CELLS, initial_partition, bias_model,
)

n_clones = size(count_matrix, 2)
@info "simulated count matrix" size = size(count_matrix) total_cells = sum(count_matrix)

# ---------------------------------------------------------------------------
# 4. Synthesise CDR3 strings, one per clonotype
# ---------------------------------------------------------------------------
const AA = collect("ACDEFGHIKLMNPQRSTVWY")

function random_cdr3(rng = Random.GLOBAL_RNG; minlen = 11, maxlen = 17)
    L = rand(rng, minlen:maxlen)
    body = String(rand(rng, AA, L - 2))
    return "C" * body * "F"
end

trb_per_clone = [random_cdr3() for _ in 1:n_clones]
tra_per_clone = [rand() < 0.55 ? random_cdr3() : "-" for _ in 1:n_clones]

# ---------------------------------------------------------------------------
# 5. Emit a long-form (one row per cell) clone table + UMAP coordinates
# ---------------------------------------------------------------------------
# Place cell-type centroids around a circle so the simulated UMAP is readable.
angles = range(0, 2π; length = N_CELL_TYPES + 1)[1:end-1]
centroids = Dict(CELL_TYPES[i] => (5cos(angles[i]), 5sin(angles[i]))
                 for i in 1:N_CELL_TYPES)

clone_rows = NamedTuple[]
umap_rows  = NamedTuple[]
bc = 0

for ci in 1:N_CELL_TYPES, ki in 1:n_clones
    n_cells = count_matrix[ci, ki]
    n_cells == 0 && continue
    cell_type = cluster_names[ci]
    clono     = "clone$(ki)"
    trb = trb_per_clone[ki]
    tra = tra_per_clone[ki]
    cdr3s = tra == "-" ? "TRB:$trb" : "TRA:$tra;TRB:$trb"
    cx, cy = centroids[cell_type]
    for _ in 1:n_cells
        bc += 1
        barcode = "BC" * lpad(bc, 6, '0') * "-1"
        push!(clone_rows, (
            Barcode = barcode,
            cell_types = cell_type,
            Clonotype = clono,
            cdr3s_aa_TRA_TRB = cdr3s,
            TRA_cdr3aa = tra,
            TRB_cdr3aa = trb,
        ))
        push!(umap_rows, (
            Barcode = barcode,
            cell_types = cell_type,
            UMAP_1 = cx + 0.7 * randn(),
            UMAP_2 = cy + 0.7 * randn(),
        ))
    end
end

# Random shuffle so the file resembles a real per-cell table
perm = Random.shuffle(1:length(clone_rows))
clone_df = DataFrame(clone_rows[perm])
umap_df  = DataFrame(umap_rows[perm])

# ---------------------------------------------------------------------------
# 6. Write to disk
# ---------------------------------------------------------------------------
out_dir = joinpath(@__DIR__, "data")
isdir(out_dir) || mkpath(out_dir)

CSV.write(joinpath(out_dir, "simulated_clone_data.tsv"),  clone_df; delim = '\t')
CSV.write(joinpath(out_dir, "simulated_UMAP_coords.tsv"), umap_df;  delim = '\t')

println("Wrote ", joinpath(out_dir, "simulated_clone_data.tsv"))
println("Wrote ", joinpath(out_dir, "simulated_UMAP_coords.tsv"))
