function ExtractNodeCoords(tree, marg_dict)
    marg_coord = Dict()
    for k in collect(keys(marg_dict))
        if !isleafnode(k)
            marg_coord[k.nodeindex] = [marg_dict[k][1].mean, marg_dict[k][2].mean];
        end
    end
    for l in getleaflist(tree)
        marg_coord[l.nodeindex] = [l.message[1].mean, l.message[2].mean];
    end
    return marg_coord
end


function PlotNode(p, node, coord)
    c = :black
    s = 4
    if !isleafnode(node)
            c = :blue
            s = 3
    end
    if MolecularEvolution.isroot(node)
            c = :red
            s = 3
    end
    scatter!(p, [coord[node.nodeindex][1]], [coord[node.nodeindex][2]], markersize=s, label="", markercolor=c,
                markerstrokewidth=0,
                marker=:circle, 
#             series_annotations=Plots.text.([node.name], :left, :bottom, :black, pointsize=10)
             )
end


"""
    PlotTreeOnUmap(tree, model, p) -> (marg_coord, anim)

Project `tree` onto the existing 2-D plot `p` (e.g. a UMAP scatter)
using the marginal node positions implied by the 2-component
`model::Vector{<:BranchModel}`. Returns the per-node marginal
coordinates and a `Plots.@animate` object that incrementally adds nodes
and branches — useful with `Plots.gif` for an animated tree walk.

See also [`PlotTreeOnUmapNoAnim`](@ref) (solid static rendering) and
[`PlotTreeOnUmapNoAnimShadow`](@ref) (faint shadow lines, for posterior
clouds).
"""
function PlotTreeOnUmap(tree, model, p)
    for (i, n) in enumerate(getnodelist(tree))
        n.nodeindex = i
    end

    marg_dict = marginal_state_dict(tree, model);
    marg_coord = ExtractNodeCoords(tree, marg_dict);

    anim = @animate for node in getnodelist(tree)
        PlotNode(p, node, marg_coord)
        for child in node.children
            if child.branchlength < 0.05
                arr = false
            else
                arr = (:closed, 0.8)
            end

            plot!(p,
                  [marg_coord[node.nodeindex][1], marg_coord[child.nodeindex][1]],
                  [marg_coord[node.nodeindex][2], marg_coord[child.nodeindex][2]],
                  markersize=8,
                  arrow=arr,
                  linewidth=2,
                  linecolor=:black, label="")
            PlotNode(p, child, marg_coord)
        end
    end

    return marg_coord, anim
end


"""
    PlotTreeOnUmapNoAnimShadow(tree, model, p) -> (marg_coord, p)

Static, transparent variant of [`PlotTreeOnUmap`](@ref). Adds the tree's
branches to the existing plot `p` as faint, alpha-blended grey lines —
intended for stacking many posterior trees on top of a UMAP scatter.
"""
function PlotTreeOnUmapNoAnimShadow(tree, model, p; linealpha = 0.1, linewidth = 0.15)
    for (i, n) in enumerate(getnodelist(tree))
        n.nodeindex = i
    end

    marg_dict = marginal_state_dict(tree, model);
    marg_coord = ExtractNodeCoords(tree, marg_dict);

    for node in getnodelist(tree)
        for child in node.children
            arr = false

            plot!(p,
                  [marg_coord[node.nodeindex][1], marg_coord[child.nodeindex][1]],
                  [marg_coord[node.nodeindex][2], marg_coord[child.nodeindex][2]],
                  marker=false,
                  arrow=arr, alpha=linealpha,
                  linewidth = linewidth,
                  linecolor="#494848", 
                  label="")
        end
    end

    return marg_coord, p
end


"""
    PlotTreeOnUmapNoAnim(tree, model, p) -> (marg_coord, p)

Static, opaque variant of [`PlotTreeOnUmap`](@ref). Draws the tree as
solid black branches with coloured nodes (root in red, internal in blue,
leaves in black) on top of the supplied 2-D plot `p` — typically used to
overlay a HIPSTR consensus tree on a UMAP.
"""
function PlotTreeOnUmapNoAnim(tree, model, p)

    for (i, n) in enumerate(getnodelist(tree))
        n.nodeindex = i
    end

    marg_dict = marginal_state_dict(tree, model);
    marg_coord = ExtractNodeCoords(tree, marg_dict);

    for node in getnodelist(tree)
        for child in node.children
            if child.branchlength < 0.05
                arr = false
            else
                arr = (:closed, 0.8)
            end

            plot!(p,
                  [marg_coord[node.nodeindex][1], marg_coord[child.nodeindex][1]],
                  [marg_coord[node.nodeindex][2], marg_coord[child.nodeindex][2]],
                  arrow=arr,
                  linewidth=2,
                  linecolor=:black, label="")
        end
    end
    for node in getnodelist(tree)
        PlotNode(p, node, marg_coord)
        for child in node.children
            PlotNode(p, child, marg_coord)
        end
    end
    return marg_coord, p
end