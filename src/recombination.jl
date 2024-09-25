"""
    function recombine(cluster_names::Vector{String}, count_matrix::Matrix{Int64}, dest::String, src::String)

Recombines cell clusters `dest` and `src`. Returns recombined cluster names and count matrix. Here recombinination means taking a clonotype-wise sum of counts from `count_matrix`. `src` and its column in the count matrix are discarded in the recombined cluster names and count matrix respectively. `dest *= "+" * src` in the recombined cluster names.

# Example
```jldoctest
julia> cluster_names = ["Type\$i" for i = 1:3];

julia> count_matrix = [i * j for i = 1:4, j = 1:3]
4×3 Matrix{Int64}:
 1  2   3
 2  4   6
 3  6   9
 4  8  12

julia> recombined_cluster_names, recbomined_count_matrix = recombine(cluster_names, count_matrix, "Type1", "Type3");

julia> recombined_cluster_names
2-element Vector{String}:
 "Type1+Type3"
 "Type2"

julia> recbomined_count_matrix
4×2 Matrix{Int64}:
  4  2
  8  4
 12  6
 16  8
```
"""
function recombine(
    cluster_names::Vector{String},
    count_matrix::Matrix{Int64},
    dest::String,
    src::String,
)
    @assert dest != src
    dest_index = findfirst(x -> x == dest, cluster_names)
    src_index = findfirst(x -> x == src, cluster_names)
    @assert !isnothing(dest_index)
    @assert !isnothing(src_index)

    recombined_cluster_names = cluster_names[Not(src_index)]
    recombined_count_matrix = count_matrix[:, Not(src_index)]
    if src_index < dest_index
        dest_index -= 1
    end
    recombined_cluster_names[dest_index] *= "+" * cluster_names[src_index]
    @views recombined_count_matrix[:, dest_index] .+= count_matrix[:, src_index]
    return recombined_cluster_names, recombined_count_matrix
end
