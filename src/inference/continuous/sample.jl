"""
    FrequencySampler

A type that allows you to specify a additive proposal function. It also holds the acceptance ratio acc_ratio (acc_ratio[1] stores the number of accepts, and acc_ratio[2] stores the number of rejects).
"""
struct FrequencySampler <: MolecularEvolution.UnivariateSampler
    acc_ratio::Array{Int64,1}
    proposal::Distribution
    FrequencySampler(proposal) = new(zeros(Int64, 2), proposal)
end

#Accept/reject individual samples
#log_posterior: Array{Float64,1} -> Array{Float64,1}
function frequencies_metropolis(
    log_posterior,
    modifier::FrequencySampler,
    curr_values::Array{Float64,1},
)
    # Adding additive normal symmetrical noise to ensure the proposal function is symmetric.
    n = length(curr_values)
    proposals = curr_values .+ rand(modifier.proposal, n)
    # The standard Metropolis acceptance criterion.
    U = rand(n)
    post_quotients = exp.(copy(log_posterior(proposals)) .- log_posterior(curr_values))
    return ifelse.(U .<= post_quotients, proposals, curr_values)
end
