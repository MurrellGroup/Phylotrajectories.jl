using Phylotrajectories, Random
using Test

@testset "Phylotrajectories.jl" begin
    @testset "importing" begin
        _, cluster_names1, count_matrix1 = import_count_matrix(
            "data/clone_data_HDM.tsv",
            :Clonotype,
            :cell_types,
            :TRB_cdr3aa,
            cluster_filters = ["Proliferating"],
        )
        _, cluster_names2, count_matrix2 = import_count_matrix("data/Clone_counts_HDM.csv")
        @test cluster_names1 == cluster_names2
        @test count_matrix1 == count_matrix2
    end

    @testset "inference" begin
        _, cluster_names, count_matrix = import_count_matrix("data/Clone_counts_HDM.csv")
        Random.seed!(1234)
        newtree1, model1, states1, LL1 = tree_inference(
            cluster_names,
            count_matrix,
            jump = 0.1,
            a = 1.0,
            b = 1.0,
            Ne = 1.0,
            rate = 50.0,
            start_branch_length = 0.1,
            max_cycles = 10,
        )
        @test LL1 ≈ -8785.964192745085
    end
end
