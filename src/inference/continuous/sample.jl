"""
    FrequencySampler(proposal::Distribution{Univariate,Continuous})

A type that allows you to specify an additive proposal function. It also holds the acceptance ratio acc_ratio (acc_ratio[1] stores the number of accepts, and acc_ratio[2] stores the number of rejects).
"""
struct FrequencySampler
    acc_ratio::Array{Int64,1}
    proposal::Distribution{Univariate,Continuous}
    FrequencySampler(proposal) = new(zeros(Int64, 2), proposal)
end

MolecularEvolution.proposal(modifier::FrequencySampler, curr_value::Array{Float64,1}) =
    curr_value + rand(modifier.proposal, length(curr_value))
MolecularEvolution.log_prior(modifier::FrequencySampler, x) = 0.0 #Constant, since, prior comes from root
function MolecularEvolution.apply_decision(modifier::FrequencySampler, accepts::BitArray)
    for b in accepts
        MolecularEvolution.apply_decision(modifier, b)
    end
end
#This should idiomatically be achieved by a loop or broadcasting metropolis_step over multiple functions and Float64s
#but the driving factor for why I specialize metropolis_step is
#it's convenient for LL to return Array{Float64,1}
function MolecularEvolution.metropolis_step(
    LL::Function,
    modifier::FrequencySampler,
    curr_value::Array{Float64,1},
)
    prop = MolecularEvolution.proposal(modifier, curr_value)
    accept_proposal =
        rand(length(curr_value)) .<=
        exp.(
            LL(prop) .+ MolecularEvolution.log_prior(modifier, prop) .- LL(curr_value) .-
            MolecularEvolution.log_prior(modifier, curr_value),
        )
    MolecularEvolution.apply_decision(modifier, accept_proposal)
    return ifelse.(accept_proposal, prop, curr_value)
end

"""
    GaussianSampler(prior::Distribution{Multivariate,Continuous}, proposal::Distribution{Multivariate,Continuous})

A type that allows you to specify multivariate proposal and prior functions over the Gaussian parameters [μ, ν] ~ N(μ, exp(ν)). It also holds the acceptance ratio `acc_ratio` (`acc_ratio[1]` stores the number of accepts, and `acc_ratio[2]` stores the number of rejects).
"""
struct GaussianSampler
    acc_ratio::Array{Int64,1}
    proposal::Distribution{Multivariate,Continuous}
    prior::Distribution{Multivariate,Continuous}
    function GaussianSampler(proposal, prior)
        @assert length(proposal) == length(prior) == 2 "Proposal and prior must have exactly 2 dimensions"
        new(zeros(Int64, 2), proposal, prior)
    end
end

# Transform mean and variance to mean and log-variance 
function tr_gaussian_params(curr_values)
    curr_values .|> [identity, log]
end

# Inverse transform mean and log-variance to mean and variance
function invtr_gaussian_params(tr_curr_values)
    tr_curr_values .|> [identity, exp]
end

#Assuming that curr_values = [μ, exp(ν)]
function MolecularEvolution.proposal(modifier::GaussianSampler, curr_values)
    return invtr_gaussian_params(tr_gaussian_params(curr_values) .+ rand(modifier.proposal))
end
function MolecularEvolution.log_prior(modifier::GaussianSampler, curr_values)
    return logpdf(modifier.prior, tr_gaussian_params(curr_values))
end
