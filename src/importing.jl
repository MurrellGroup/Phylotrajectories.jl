function create_dict_from_dataframe(df::DataFrame, key_col::Symbol, value_col::Symbol)
    dict = Dict{String,Vector{String}}()
    for row in eachrow(df)
        key = row[key_col]
        value = row[value_col]
        if haskey(dict, key)
            push!(dict[key], value)
        else
            dict[key] = [value]
        end
    end
    return dict
end

function import_count_dataframe(
    data::DataFrame,
    clone_column_name::Symbol,
    clusters_column_name::Symbol,
    cdr3_column_name::Symbol,
)
    # Force these to be strings for invariance purposes
    transform!(
        data,
        [clone_column_name, clusters_column_name, cdr3_column_name] .=> ByRow(string),
        renamecols = false, #Possible without this?
    )
    # Count the occurrences of each combination of cell_type and clonotype
    count_df =
        combine(groupby(data, [clusters_column_name, clone_column_name]), nrow => :count)
    # Pivot the DataFrame to get the desired format
    pivot_df = unstack(count_df, clone_column_name, :count)
    # Replace missing values with 0
    pivot_df = coalesce.(pivot_df, 0)
    # Transposing and removing the first column that includes clones names
    pivot_df = permutedims(pivot_df, 1)

    cdr3_dict = create_dict_from_dataframe(
        unique(data[:, [clone_column_name, cdr3_column_name]]),
        clone_column_name,
        cdr3_column_name,
    )
    pivot_df[!, cdr3_column_name] =
        map(x -> cdr3_dict[x], pivot_df[:, clusters_column_name])

    return pivot_df[:, [clusters_column_name, cdr3_column_name]],
    pivot_df[:, Not([clusters_column_name, cdr3_column_name])]
end

function import_count_matrix(
    fname,
    clone_column_name::Symbol,
    clusters_column_name::Symbol,
    cdr3_column_name::Symbol;
    cluster_filters::Vector{String} = [""],
)
    df = CSV.read(fname, delim = '\t', DataFrame)
    df = df[df[!, clusters_column_name].∉Ref(cluster_filters), :]
    clono_info, count_matrix_df = import_count_dataframe(
        df,
        clone_column_name,
        clusters_column_name,
        cdr3_column_name,
    )
    cluster_sizes = Dict(eachrow(combine(groupby(df, [clusters_column_name]), nrow => :count)))
    return clono_info, names(count_matrix_df), cluster_sizes, Matrix(Matrix(count_matrix_df)'), df #Second Matrix call needed?
end

function import_count_matrix(fname)
    count_matrix_df = CSV.read(fname, DataFrame)
    return nothing, names(count_matrix_df), Matrix(Matrix(count_matrix_df)') #Second Matrix call needed?
end
