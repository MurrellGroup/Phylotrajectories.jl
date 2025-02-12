#TODO: consider slipping in more @assert for dimension matching

#-----------------Gaussian prop---------------
#Partition behavior
mutable struct IndependentGaussiansPartition <: MolecularEvolution.ContinuousPartition
    means::Array{Float64,1}
    vars::Array{Float64,1}
    norm_consts::Array{Float64,1}
    counts::Array{Int64,1}
    function IndependentGaussiansPartition(means, vars, norm_consts)
        @assert length(means) == length(vars) == length(norm_consts)
        new(means, vars, norm_consts)
    end
end

IndependentGaussiansPartition(means, vars) = new(means, vars, zeros(length(means)))

IndependentGaussiansPartition(n) =
    IndependentGaussiansPartition(zeros(n), ones(n), zeros(n))

Base.length(g::IndependentGaussiansPartition) = length(g.means)

function Base.getindex(g::IndependentGaussiansPartition, i::Int)
    1 <= i <= length(g) || throw(BoundsError(g, i))
    return g.means[i], g.vars[i], g.norm_consts[i]
end
function Base.setindex!(
    g::IndependentGaussiansPartition,
    v::Tuple{Float64,Float64,Float64},
    i::Int,
)
    1 <= i <= length(g) || throw(BoundsError(g, i))
    g.means[i] = v[1]
    g.vars[i] = v[2]
    g.norm_consts[i] = v[3]
end

#Overloading the copy_partition to avoid deepcopy.
function MolecularEvolution.copy_partition(src::IndependentGaussiansPartition)
    return IndependentGaussiansPartition(
        copy(src.means),
        copy(src.vars),
        copy(src.norm_consts),
    )
end

function MolecularEvolution.copy_partition_to!(
    dest::IndependentGaussiansPartition,
    src::IndependentGaussiansPartition,
)
    dest.means .= src.means
    dest.vars .= src.vars
    dest.norm_consts .= src.norm_consts
end

#=
For merge_two_gaussians and gaussian_pdf, consider checking for edge cases first (as in GaussianPartition),
use a filter for the computations
=#


#From the first section of http://www.tina-vision.net/docs/memos/2003-003.pdf
function MolecularEvolution.merge_two_gaussians(
    g1::IndependentGaussiansPartition,
    g2::IndependentGaussiansPartition,
)
    res_gaussians = IndependentGaussiansPartition(length(g1))
    res_gaussians.vars .= 1 ./ (1 ./ g1.vars .+ 1 ./ g2.vars)
    res_gaussians.means .=
        res_gaussians.vars .* (g1.means ./ g1.vars .+ g2.means ./ g2.vars)
    # log of scaling constant
    res_gaussians.norm_consts .=
        -0.5 .* (
            log.(2 .* pi * (g1.vars .* g2.vars ./ res_gaussians.vars)) .+
            (g1.means .^ 2 ./ g1.vars) .+ (g2.means .^ 2 ./ g2.vars) .-
            (res_gaussians.means .^ 2 ./ res_gaussians.vars)
        )
    res_gaussians.norm_consts .+= (g1.norm_consts .+ g2.norm_consts)
    for i = 1:length(g1)
        #Handling some edge cases. These aren't mathematically sensible. A gaussian with "Inf" variance will behave like a 1,1,1,1 vector in discrete felsenstein.
        #To-do: update some of these so that the norm constant is properly handled, even if the variance is Inf (so it isn't exactly well-defined anyway)
        if g1.vars[i] == 0 && g2.vars[i] == 0 && g1.means[i] != g2.means[i]
            error("both gaussians have 0 variance but different means")
        elseif g1.vars[i] == 0
            res_gaussians[i] = MolecularEvolution._merge_point_mass(g1, g2, i)
        elseif g2.vars[i] == 0
            res_gaussians[i] = MolecularEvolution._merge_point_mass(g2, g1, i)
        end
        if g1.vars[i] == Inf && g2.vars[i] == Inf
            res_gaussians[i] = ((g1.means[i] + g2.means[i]) / 2, Inf, 0.0)
        elseif g1.vars[i] == Inf
            res_gaussians[i] = g2[i]
        elseif g2.vars[i] == Inf
            res_gaussians[i] = g1[i]
        end
    end
    return res_gaussians
end

function MolecularEvolution._merge_point_mass(point::IndependentGaussiansPartition, regular::IndependentGaussiansPartition, i::Int)
    mean, var, norm_const = point[i]
    norm_const += logpdf(Normal(mean, sqrt(var)), regular.means[i]) + regular.norm_consts[i]
    return mean, var, norm_const
end

function MolecularEvolution.identity!(dest::IndependentGaussiansPartition)
    dest.vars .= Inf
    dest.norm_consts .= 0.0
end

function MolecularEvolution.combine!(
    dest::IndependentGaussiansPartition,
    src::IndependentGaussiansPartition,
)
    new_g = MolecularEvolution.merge_two_gaussians(dest, src)
    dest.means .= new_g.means
    dest.vars .= new_g.vars
    dest.norm_consts .= new_g.norm_consts
end

function MolecularEvolution.site_LLs(part::IndependentGaussiansPartition)
    return part.norm_consts
end

function MolecularEvolution.gaussian_pdf(
    g::IndependentGaussiansPartition,
    x::Array{Float64,1},
)
    result = pdf.(Normal.(g.means, sqrt.(g.vars)), x)
    for i = 1:length(g)
        if g.vars[i] == 0
            result[i] = Float64(x == g.means[i]) #Hokey...
        end
    end
    return result
end

#And sampling
function MolecularEvolution.sample_partition!(partition::IndependentGaussiansPartition)
    partition.means .= randn(length(partition)) .* sqrt.(partition.vars) .+ partition.means
    partition.vars .= 0.0
    partition.norm_consts .= 0.0
end

#And max
function MolecularEvolution.max_partition!(partition::IndependentGaussiansPartition)
    partition.vars .= 0.0
    partition.norm_consts .= 0.0
end

function MolecularEvolution.obs2partition!(
    partition::IndependentGaussiansPartition,
    obs::Array{Float64,1},
)
    partition.means .= obs
    partition.vars .= 0.0
    partition.norm_consts .= 0.0
end

function MolecularEvolution.partition2obs(partition::IndependentGaussiansPartition)
    return partition.means
end
