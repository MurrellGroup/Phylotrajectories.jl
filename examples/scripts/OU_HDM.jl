import Pkg
Pkg.activate("/home/masha/Phylotrajectories.jl")

using Phylotrajectories
using MolecularEvolution

using CSV, DataFrames, StatsBase, Distributions, Phylo, Plots, Statistics, LinearAlgebra, Random, ForwardBackward

using ProgressMeter
using SpecialFunctions
using Plots.PlotMeasures

using Plots
using ColorSchemes
using Printf

using Random
number_seed = 1234
Random.seed!(number_seed)

path = "/home/masha/CLEAN_PhyloTrajectory_paper/results/HDM"

clono_info, cluster_names, _,count_matrix = import_count_matrix(
    "/home/masha/PhyloTraj_paper/Final_notebooks/HDM/HDM_clone_data_l2_notcr_091125.tsv",
    :Clonotype,
    :cell_types_l2,
    :TRB_cdr3aa,
    cluster_filters = ["Proliferating"] 
);


################### Inference ####################

digamma_value = 0.5
trigamma_value = 0.5

eqmu = 1.5
eqtheta = 0.1
v = 1.0

tree_warmup_cycles = 100
start_branch_length = 0.1

function make_ou_model(;
    nni = 1,
    branchlength = 1,
    root = 1,
    models = 1,
    tree_warmup_cycles = 100,
    burn_in = 10000,
    sample_interval = 100,
    n_samples = 100,
)
    return OUContinuousModel(
        update = OUContinuousUpdate(
            nni = nni,
            branchlength = branchlength,
            root = root,
            models = models,
            branchlength_sampler = BranchlengthSampler(Normal(0, 0.1), Normal(-1, 1)),
            root_sampler = OUGaussianStateSample(MvNormal(zeros(2), Diagonal([0.01, 0.01])), 
                            MvNormal(zeros(2), Diagonal([1.0, 0.1])), 
                1e-1, 1),
            ou_eqmu_sampler = OUEqmuSampler(Normal(0.0, 2.0), Normal(1.5, 1.0), 1.0, 0.1),
            ou_theta_sampler = OUThetaSampler(Normal(0, 0.5), Normal(-1, 1), 1.5, 1.0),
        ),
        start_branch_length = start_branch_length,
        tree_warmup_cycles = tree_warmup_cycles,
        burn_in = burn_in,
        sample_interval = sample_interval,
        n_samples = n_samples,
    )
end

println("="^60)
println("Running OU MCMC inference...")
println("="^60)

plot_init, init_tree, trees, LLs, models, root_ps, upd =
    tree_inference(
        make_ou_model(
            burn_in = 10000,
            sample_interval = 100,
            n_samples = 1000,
            tree_warmup_cycles = 100
        ),
    cluster_names, count_matrix;
    eqmu=eqmu, eqtheta=eqtheta, v=v,
    d=digamma_value, g=trigamma_value,
)

####################### Plots and saving results #######################

ladderize!.(trees)
hip, node2logcred, node2support = HIPSTR(trees; getcred=true, getsupport=true)
mltpl = plot_multiple_trees(trees, hip; line_width=0.075);

lls = plot(LLs, legend=:none);
title!(lls, "LLs");

means = plot([x[3] for x in models], legend=:none);
title!(means, ("EqMu Acc Ratio $(round(upd.models_update.eqmu_sampler.acc_ratio[1], digits=3))"));

theta = plot([x[1] for x in models], legend=:none);
title!(theta, ("Theta Acc Ratio $(round(upd.models_update.theta_sampler.acc_ratio[1], digits=3))"));

vars = plot([x[2] for x in models], legend=:none);
title!(vars, ("Var"));

plot_root = plot(xlims=(-5, 5), legend=false)
for p in root_ps
    plot!(plot_root, x->pdf(Normal(p[1][1], p[1][2]), x))
end
title!(plot_root, "Root Acc Position: $(round(upd.root_update.acc_ratio.position[1], digits=3)) State: $(round(upd.root_update.acc_ratio.state[1], digits=3))");

combined_plt = plot(
    [means, vars, theta, lls, plot_root, mltpl]...,
    layout = (2, 3), size = (1500, 600)
    );

hipster_tree = plot(mltpl, fontsize=10, margins=10Plots.PlotMeasures.mm);

savefig(combined_plt, "$(path)/OU_HDM_$(digamma_value)_$(trigamma_value)_$(number_seed)_theta.pdf")
savefig(hipster_tree, "$(path)/OU_HDM_$(digamma_value)_$(trigamma_value)_hipster_$(number_seed)_theta.pdf")

# Save a representative tree from the MCMC samples (e.g., the last tree)
Phylo.write("$(path)/OU_HDM_$(digamma_value)_$(trigamma_value)_HIPSTR_$(number_seed)_theta.newick", get_phylo_tree(hip))

####################### Clone clustering #######################


function max_root_distance(node, dist = 0.0)
    if isleafnode(node)
        return dist
    end
    return maximum(max_root_distance(ch, dist + ch.branchlength) for ch in node.children)
end


for n in getnodelist(hip)
    if isleafnode(n)
        n.node_data = Dict(
            "size" => 2.5,
            "support_plot" => 0.0,
        )
    else
        s = node2support[n]
        n.node_data = Dict(
            "size" => 8 + 18 * sqrt(s),
            "support_plot" => s,
        )
    end
end

ladderize!(hip)
phylo_tree = get_phylo_tree(hip)

max_depth = max_root_distance(hip)
xpad_left = 0.06 * max_depth
xpad_right = 0.28 * max_depth

pl_cred = plot(
    phylo_tree,
    showtips = true,
    tipfont = 10,
    markersize = values_from_phylo_tree(phylo_tree, "size"),
    marker_z = values_from_phylo_tree(phylo_tree, "support_plot"),
    markercolor = cgrad(reverse(ColorSchemes.RdBu_10.colors)),
    clims = (0.0, 1.0),
    markerstrokecolor = :black,
    markerstrokewidth = 0.25,
    linewidth = 1.0,
    xlims = (-xpad_left, max_depth + xpad_right),
    left_margin = 10mm,
    right_margin = 35mm,
    top_margin = 8mm,
    bottom_margin = 8mm,
    colorbar = :right,
    size = (950, 700),
)

for (i, s) in enumerate(pl_cred.series_list)
    println("Series $i: type=$(s.plotattributes[:seriestype]), npoints=$(length(s[:x]))")
end

scatter_series = pl_cred.series_list[2]
xs = scatter_series[:x]
ys = scatter_series[:y]

zs = scatter_series[:marker_z]

for i in eachindex(xs)
    s = zs[i]
    if s > 0.0  # skip leaves (which you set to 0.0)
        label = string(round(s; digits = 3))
        annotate!(pl_cred,
            xs[i] + 0.07 * max_depth,  # offset right of the bubble
            ys[i],
            text(label, 8, :left, :black)
        )
    end
end

savefig(pl_cred, "$(path)/OU_HDM_$(digamma_value)_$(trigamma_value)_HIPSTR_with_CRED_seed_$(number_seed).pdf");


run_ou_and_build_clone_matrix(trees,
    models,
    count_matrix,
    cluster_names,
    0.5,
    0.5,
    "$(path)/OU_HDM_0.5_0.5_clone_matrix_seed_$(number_seed).csv"
    )

UMAP_coords = CSV.read("/home/masha/PhyloTraj_paper/Final_notebooks/HDM/UMAP_coords_l2_091125.tsv", delim='\t',DataFrame);
UMAP_coords.UMAP_1 = UMAP_coords.UMAP_1 .* (-1);

UMAP_coords_gdf = groupby(UMAP_coords, :cell_types_l2);
centroids_umap = combine(UMAP_coords_gdf, [:UMAP_1, :UMAP_2] => ((u1, u2) -> (C1=mean(u1), C2=mean(u2))) => AsTable);

color_palette = Dict(
    "Th2_2"         => "#1a7fe4",  # Blue
    "Naive"         => "#49ce50",  # Green
    "Act_circ"      => "#E5D156",  # Yellow
    "Act(resid)"    => "#8A4AFA",  # Bright violet
    "Tregs-0"       => "#E74C3C",  # True red
    "Tregs-1"       => "#00BCD4",  # Cyan
    "Tregs-2"       => "#207541",  # Orange
    "Tregs-3"       => "#A4008B",  # Pink
    "Th2"           => "#708090",  # Slate grey
    "IFN_response"  => "#8C564B",  # Brown
    "Proliferating" => "#FFBDF4"   # Red-orange
)

sc_p = scatter(UMAP_coords.UMAP_1, 
        UMAP_coords.UMAP_2,
        group=UMAP_coords.cell_types_l2,
        color=[color_palette[c] for c in UMAP_coords.cell_types_l2], 
        markerstrokewidth=0, 
        alpha=0.4);

for tree in trees
    newtree_umap = copy_tree(tree);
    message_template_umap = [GaussianPartition(0., 1000.), GaussianPartition(0., 1000.)];

    bm_model = [MolecularEvolution.BrownianMotion(0.,0.5), MolecularEvolution.BrownianMotion(0.,0.5)];
    internal_message_init!(newtree_umap, message_template_umap);

    for x in getleaflist(newtree_umap)
        tmp = filter(row->row[:cell_types_l2] == x.name, centroids_umap);
        x.message = [GaussianPartition(tmp.C1[1], 0.0), GaussianPartition(tmp.C2[1], 0.0)];
    end

    felsenstein!(newtree_umap, bm_model);
    global sc_p
    _, sc_p = PlotTreeOnUmapNoAnimShadow(newtree_umap, bm_model, sc_p);
end


newtree_umap = copy_tree(hip);
message_template_umap = [GaussianPartition(0., 1000.), GaussianPartition(0., 1000.)]

bm_model = [MolecularEvolution.BrownianMotion(0.,0.5), MolecularEvolution.BrownianMotion(0.,0.5)]
internal_message_init!(newtree_umap, message_template_umap) 

for x in getleaflist(newtree_umap)
    tmp = filter(row->row[:cell_types_l2] == x.name, centroids_umap)
    x.message = [GaussianPartition(tmp.C1[1], 0.0), GaussianPartition(tmp.C2[1], 0.0)]
end

felsenstein!(newtree_umap, bm_model)

coord, im_with_hip = PlotTreeOnUmapNoAnim(newtree_umap, bm_model, sc_p);

savefig(im_with_hip, "$(path)/HDM_with_MCMC_trees.pdf")

sc_p = scatter(UMAP_coords.UMAP_1, 
        UMAP_coords.UMAP_2,
        group=UMAP_coords.cell_types_l2,
        color=[color_palette[c] for c in UMAP_coords.cell_types_l2], 
        markerstrokewidth=0, 
        alpha=0.4);

newtree_umap = copy_tree(hip);
message_template_umap = [GaussianPartition(0., 1000.), GaussianPartition(0., 1000.)]

bm_model = [MolecularEvolution.BrownianMotion(0.,0.5), MolecularEvolution.BrownianMotion(0.,0.5)]
internal_message_init!(newtree_umap, message_template_umap) 

for x in getleaflist(newtree_umap)
    tmp = filter(row->row[:cell_types_l2] == x.name, centroids_umap)
    x.message = [GaussianPartition(tmp.C1[1], 0.0), GaussianPartition(tmp.C2[1], 0.0)]
end

felsenstein!(newtree_umap, bm_model)

coord, im_with_hip = PlotTreeOnUmapNoAnim(newtree_umap, bm_model, sc_p);

savefig(im_with_hip, "$(path)/HDM_no_MCMC_trees.pdf")


sc_p = scatter(UMAP_coords.UMAP_1, 
        UMAP_coords.UMAP_2,
        group=UMAP_coords.cell_types_l2,
        color=[color_palette[c] for c in UMAP_coords.cell_types_l2], 
        markerstrokewidth=0, 
        alpha=0.4);

newtree_umap = copy_tree(hip);
message_template_umap = [GaussianPartition(0., 1000.), GaussianPartition(0., 1000.)]

bm_model = [MolecularEvolution.BrownianMotion(0.,0.5), MolecularEvolution.BrownianMotion(0.,0.5)]
internal_message_init!(newtree_umap, message_template_umap) 

for x in getleaflist(newtree_umap)
    tmp = filter(row->row[:cell_types_l2] == x.name, centroids_umap)
    x.message = [GaussianPartition(tmp.C1[1], 0.0), GaussianPartition(tmp.C2[1], 0.0)]
end

felsenstein!(newtree_umap, bm_model)

coord, im_with_hip_anim = PlotTreeOnUmap(newtree_umap, bm_model, sc_p);

gif(im_with_hip_anim, "$(path)/HDM_HIPSTR.gif", fps = 0.5)