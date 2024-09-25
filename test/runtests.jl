using Phylotrajectories, MolecularEvolution, Random
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
        newtree1, model1, states1, LL1, LLs1 = tree_inference(
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

        tree_inference(
            cluster_names,
            count_matrix,
            jump = 0.1,
            a = 1.0,
            b = 1.0,
            Ne = 1.0,
            rate = 50.0,
            start_branch_length = 0.1,
            max_cycles = 10,
            n_random_trees = 2,
        )
    end

    @testset "simulations" begin
        Random.seed!(4321)
        n_cell_types = 10
        n_clonotypes = 10
        n_cells = 100 * n_cell_types * n_clonotypes
        n(t) = (10 * n_cell_types) / (1 + exp(t - 10))
        tree = sim_tree(n_cell_types, n, n_cell_types / 5, mutation_rate = 0.05)
        initial_partition = GaussianPartition(0.0, 0.0)
        default_BM_model = BrownianMotion(-0.3, 1.5)

        #Bias model
        bias_clonotype = 3
        towards_cluster = 1
        pos_bias_BM_model = BrownianMotion(1.5, 1.5)
        neg_bias_BM_model = BrownianMotion(-10.0, 1.5)
        function bias_model(i)
            if i != bias_clonotype
                return default_BM_model
            end
            bias_branches =
                Set(Phylotrajectories.getnode2rootpath(getleaflist(tree)[towards_cluster]))
            d = Dict{FelNode,BrownianMotion}()
            for n in getnodelist(tree)
                if n ∈ bias_branches
                    d[n] = pos_bias_BM_model
                else
                    d[n] = neg_bias_BM_model
                end
            end
            return n::FelNode -> [d[n]]
        end

        cluster_names, count_matrix = sim_count_matrix(
            tree,
            n_clonotypes,
            n_cells,
            initial_partition,
            default_BM_model,
        )
        cluster_names, count_matrix = sim_count_matrix(
            tree,
            n_clonotypes,
            n_cells,
            initial_partition,
            [BrownianMotion(-0.3, i / 2) for i = 1:n_clonotypes],
        )
        cluster_names, count_matrix =
            sim_count_matrix(tree, n_clonotypes, n_cells, initial_partition, bias_model)
    end

    @testset "recombination" begin
        cluster_names = ["Type$i" for i = 1:3]
        count_matrix = [i * j for i = 1:4, j = 1:3]
        recombined_cluster_names, recbomined_count_matrix =
            recombine(cluster_names, count_matrix, "Type1", "Type3")
        @test recombined_cluster_names == ["Type1+Type3", "Type2"]
        @test recbomined_count_matrix == [4 2; 8 4; 12 6; 16 8]
    end
end
