using Phylotrajectories, MolecularEvolution, Random, Distributions
using Test

@testset "Phylotrajectories.jl" begin
    @testset "importing" begin
        _, cluster_names1, _, count_matrix1 = import_count_matrix(
            "data/clone_data_HDM.tsv",
            :Clonotype,
            :cell_types,
            :TRB_cdr3aa,
            cluster_filters = ["Proliferating"],
        )
        _, cluster_names2, _, count_matrix2 = import_count_matrix("data/Clone_counts_HDM.csv")
        @test cluster_names1 == cluster_names2
        @test count_matrix1 == count_matrix2

        _, cluster_names3, _, count_matrix3 = import_count_matrix(
            "data/HDM_clone_data_l1_v2_1210_notcr_tmp.tsv",
            :Clonotype_tmp, #Consists of ints, want to be able to read this kind of data
            :cell_types,
            :TRB_cdr3aa,
        )
        @test count_matrix3[1:5, 1:5] ==
              [1 1 0 0 0; 0 9 0 2 0; 0 2 1 2 6; 0 3 0 1 1; 0 0 0 0 2]
    end

    @testset "inference" begin
        _, cluster_names, _, count_matrix = import_count_matrix("data/Clone_counts_HDM.csv")
        Random.seed!(1234)

        model = DiscreteModel(
            jump = 0.1,
            a = 1.0,
            b = 1.0,
            Ne = 1.0,
            sample_rate = 50.0,
            start_branch_length = 0.1,
            max_cycles = 10,
        )
        newtree1, model1, states1, trees1, LLs1 =
            tree_inference(model, cluster_names, count_matrix)
        @test maximum(LLs1) ≈ -8785.964877067909

        model_multi = DiscreteModel(
            jump = 0.1,
            a = 1.0,
            b = 1.0,
            Ne = 1.0,
            sample_rate = 50.0,
            start_branch_length = 0.1,
            max_cycles = 10,
            n_random_trees = 2,
        )
        tree_inference(model_multi, cluster_names, count_matrix)

        model_mcmc =
            DiscreteModel(ML = false, n_samples = 10, burn_in = 10, sample_interval = 10)
        tree_inference(model_mcmc, cluster_names, count_matrix)

        model_cont = ContinuousModel(n_samples = 10, burn_in = 10, sample_interval = 10)
        newtree2, trees, LLs2, models =
            tree_inference(model_cont, cluster_names, count_matrix)

        @testset "Continuous" begin
            # Init GaussianPartition and IndependentGaussiansPartition
            n = 10
            means = collect(1:n) .- 5
            vars = vcat(1:(n-1), Inf)
            norm_consts = zeros(n)
            shuffle!(means)
            idg1 = IndependentGaussiansPartition(copy(means), copy(vars), copy(norm_consts))
            gs1 = GaussianPartition.(copy(means), copy(vars), copy(norm_consts))
            shuffle!(means)
            idg2 = IndependentGaussiansPartition(copy(means), copy(vars), copy(norm_consts))
            gs2 = GaussianPartition.(copy(means), copy(vars), copy(norm_consts))

            # Test Vector{GaussianPartition} <=> IndependentGaussiansPartition
            idg3 = copy_partition(idg2)
            gs3 = copy_partition.(gs2)
            combine!(idg3, idg2)
            combine!.(gs3, gs2)
            combine!(idg1, idg2)
            combine!.(gs1, gs2)

            for i = 1:n
                @test idg1[i] == (gs1[i].mean, gs1[i].var, gs1[i].norm_const)
                @test idg3[i] == (gs3[i].mean, gs3[i].var, gs3[i].norm_const)
            end

            x = Vector{Float64}(1:n)
            @test Distributions.logpdf(idg1, x) ==
                  Distributions.logpdf.(gs1, x)
            MolecularEvolution.identity!(idg1)
            MolecularEvolution.site_LLs(idg1)
            sample_partition!(idg1)
            obs2partition!(idg1, x)
            partition2obs(idg1)

            node = FelNode(0.4, "20")
            mean_drifts = randn(n)
            var_drifts = randn(n)
            bm = BrownianMotion(0.0, 0.0)
            bms = BrownianMotion.(mean_drifts, var_drifts)
            ibm = IndependentBrownianMotion(0.0, 0.0)
            ibms = IndependentBrownianMotion(mean_drifts, var_drifts)
            forward!(idg3, idg2, ibm, node)
            backward!(idg1, idg2, ibms, node)
            for i = 1:n
                forward!(gs3[i], gs2[i], bm, node)
                backward!(gs1[i], gs2[i], bms[i], node)
                #Test Vector{GaussianPartition} + Vector{BrownianMotion} <=> IndependentGaussiansPartition + IndependentBrownianMotion
                @test idg1[i] == (gs1[i].mean, gs1[i].var, gs1[i].norm_const)
                #Test Vector{GaussianPartition} + BrownianMotion <=> IndependentGaussiansPartition + IndependentBrownianMotion (scalar)
                @test idg3[i] == (gs3[i].mean, gs3[i].var, gs3[i].norm_const)
            end
        end
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
