include("continuous/sample.jl") #we're importing this to make sure sampling types are defined

abstract type InferenceModel end

const CANONICAL_JUMP = 0.1


"""
    DiscreteModel(; <keyword arguments>)

Frequencies are discretized into states, and evolve by a continuous time Markov chain (CTMC). You may choose between searching for the Maximum Likelihood tree, or sampling trees from the posterior distribution using the Metropolis algorithm.

# Keyword Arguments
- `ML::Bool=true`: whether to perform Maximum Likelihood inference or MCMC.
- `jump::Float64=CANONICAL_JUMP`: the step size between frequency states (log domain).
- `a::Float64=1.0`: the rate of transitions to lower frequencies.
- `b::Float64=1.0`: the rate of transitions to higher frequencies.
- `Ne::Float64=1.0`: the initial tree's effective population size.
- `sample_rate::Float64=10.0`: the initial tree's sample rate.
- `start_branch_length::Float64=0.1`: the initial tree's non-root branch lengths.
- `max_cycles::Int=10`: the number of topology-only optimization iterations.
- `n_random_trees::Int=1`: the number of initial tree n_samples (for Maximum Likelihood inference).
- `n_samples::Int=10`: number of MCMC samples to collect.
- `burn_in::Int=1000`: the number of MCMC iterations to discard as burn-in.
- `sample_interval::Int=10`: the number of MCMC iterations between samples.
"""
struct DiscreteModel <: InferenceModel
    ML::Bool
    jump::Float64
    a::Float64
    b::Float64
    Ne::Float64
    sample_rate::Float64
    start_branch_length::Float64
    max_cycles::Int
    n_random_trees::Int
    n_samples::Int
    burn_in::Int
    sample_interval::Int

    function DiscreteModel(;
        ML = true,
        jump = CANONICAL_JUMP,
        a = 1.0,
        b = 1.0,
        Ne = 1.0,
        sample_rate = 10.0,
        start_branch_length = 0.1,
        max_cycles = 10,
        n_random_trees = 1,
        n_samples = 10,
        burn_in = 1000,
        sample_interval = 10,
    )
        new(ML, jump, a, b, Ne, sample_rate, start_branch_length, max_cycles, n_random_trees, n_samples, burn_in, sample_interval)
    end
end

"""
    ContinuousModel(; <keyword arguments>)

Frequencies diffuse throughout the tree in a continuous space per Brownian motion. The posterior distribution over trees is explored with the Metropolis algorithm, where the root distribution and the the frequencies at the leaves are sampled.

# Keyword Arguments
- `mean_drift::Float64=0.0`: mean drift parameter for the Brownian motion process in the unbounded frequency space.
- `Ne::Float64=1.0`: the initial tree's effective population size.
- `sample_rate::Float64=10.0`: the initial tree's sample rate.
- `start_branch_length::Float64=0.1`: the initial tree's non-root branch lengths.
- `n_samples::Int=10`: number of MCMC samples to collect after burn-in.
- `burn_in::Int=1000`: the number of MCMC iterations to discard as burn-in.
- `sample_interval::Int=10`: the number of MCMC iterations between samples.
- `consecutive_root_samples::Int=10`: number of consecutive root distribution proposals per MCMC iteration.
- `frequency_sampler::FrequencySampler=FrequencySampler(Normal())`: the proposal distribution for frequency updates in MCMC.
- `root_distribution_sampler::GaussianSampler=GaussianSampler(MvNormal(zeros(2), Diagonal([1.0, 0.1])), MvNormal(zeros(2), Diagonal([0.1, 0.1])))`: the proposal and prior distributions for root updates in MCMC.
"""
struct ContinuousModel <: InferenceModel
    mean_drift::Float64
    frequency_sampler::FrequencySampler
    root_distribution_sampler::GaussianSampler
    Ne::Float64
    sample_rate::Float64
    start_branch_length::Float64
    n_samples::Int
    burn_in::Int
    sample_interval::Int
    consecutive_root_samples::Int

    function ContinuousModel(;
        mean_drift = 0.0,
        frequency_sampler = FrequencySampler(Normal()),
        root_distribution_sampler = GaussianSampler(
            MvNormal(zeros(2), Diagonal([1.0, 0.1])),
            MvNormal(zeros(2), Diagonal([0.1, 0.1])),
        ),
        Ne = 1.0,
        sample_rate = 10.0,
        start_branch_length = 0.1,
        n_samples = 10,
        burn_in = 1000,
        sample_interval = 10,
        consecutive_root_samples = 10,
    )
        new(
            mean_drift,
            frequency_sampler,
            root_distribution_sampler,
            Ne,
            sample_rate,
            start_branch_length,
            n_samples,
            burn_in,
            sample_interval,
            consecutive_root_samples,
        )
    end
end
