##########################################################################
# Precoding types
type PrecodingResults
    results::Dict{ASCIIString, Any}
end
PrecodingResults() =
    PrecodingResults(Dict{ASCIIString, Any}())
Base.getindex(p::PrecodingResults, k::ASCIIString) =
    getindex(p.results, k)
Base.setindex!(p::PrecodingResults, v, k::ASCIIString) =
    setindex!(p.results, v, k)

# AuxPrecodingParams is defined in network.jl due to cyclic dependencies
# between files.

##########################################################################
# Reference implementation of standard coordinated precoding algorithms
include("Eigenprecoding.jl")
include("Gomadam2008_MaxSINR.jl")
include("Komulainen2013_WMMSE.jl")
include("Razaviyayn2013_MinMaxWMMSE.jl")
include("Shi2011_WMMSE.jl")
include("SumMSEMinimization.jl")

##########################################################################
# Standard functions to calculate rates from optimal MSE weights
ReferenceImplementationState = Union(
    Gomadam2008_MaxSINRState,
    Komulainen2013_WMMSEState,
    Razaviyayn2013_MinMaxWMMSEState,
    Shi2011_WMMSEState,
    SumMSEMinimizationState,
)

function calculate_logdet_rates(state::ReferenceImplementationState)
    K = length(state.W)
    ds = Int[ size(state.W[k], 1) for k = 1:K ]; max_d = maximum(ds)

    logdet_rates = zeros(Float64, K, max_d)
    for k = 1:K; if ds[k] > 0
        # W is p.d., so we should only get real eigenvalues. Numerically we may
        # get some imaginary noise however. Also, numerically the eigenvalues
        # may be less than 1, so we need to handle that to not get negative
        # rates.
        r = log2(max(1, real(eigvals(state.W[k]))))

        if ds[k] < max_d
            logdet_rates[k,:] = cat(1, r, zeros(Float64, max_d - ds[k]))
        else
            logdet_rates[k,:] = r
        end
    end; end

    return logdet_rates
end

function calculate_MMSE_rates(state::ReferenceImplementationState)
    K = length(state.W)
    ds = Int[ size(state.W[k], 1) for k = 1:K ]; max_d = maximum(ds)

    MMSE_rates = zeros(Float64, K, max_d)
    for k = 1:K; if ds[k] > 0
        # Invert W to get MMSE matrix and pick out diagonal elements, to obtain
        # MSE performance of MMSE receiver. Then take log2 of the reciprocal
        # to get the rates for the streams.
        # N.B. it is a little wasteful to take the inverse here, but it makes
        # for cleaner code in the algorithms. Also, these matrices are typically
        # small (like 2-by-2 or 3-by-3), so the complexity isn't too bad.
        E = state.W[k]\eye(state.W[k])
        r = log2(max(1, real(1./diag(E))))

        if ds[k] < max_d
            MMSE_rates[k,:] = cat(1, r, zeros(Float64, max_d - ds[k]))
        else
            MMSE_rates[k,:] = r
        end
    end; end

    return MMSE_rates
end

function calculate_allocated_power(state::ReferenceImplementationState)
    K = length(state.W)
    ds = Int[ size(state.W[k], 1) for k = 1:K ]; max_d = maximum(ds)

    allocated_power = zeros(Float64, K, max_d)

    for k = 1:K
        p = [ vecnorm(state.V[k][:,n])^2 for n = 1:ds[k] ]

        if ds[k] < max_d
            allocated_power[k,:] = cat(1, p, zeros(Float64, max_d - ds[k]))
        else
            allocated_power[k,:] = p
        end
    end

    return allocated_power
end

##########################################################################
# Initialization functions
zero_receivers(channel::SinglecarrierChannel, ds::Vector{Int}) = 
    [ zeros(Complex128, channel.Ns[k], ds[k]) for k = 1:channel.K ]

unity_MSE_weights(ds::Vector{Int}) =
    [ Hermitian(eye(ds[k])) for k = 1:length(ds) ]

function initial_precoders(channel::SinglecarrierChannel, Ps::Vector{Float64},
    sigma2s::Vector{Float64}, ds::Vector{Int}, cell_assignment::CellAssignment,
    aux_params::AuxPrecodingParams)

    V = Array(Matrix{Complex128}, channel.K)

    if aux_params["initial_precoders"] == "dft"
        for i = 1:channel.I
            served = served_MS_ids(i, cell_assignment)
            Kc = length(served)

            for k in served
                V[k] = sqrt(Ps[i]/(channel.Ms[i]*ds[k]*Kc))*fft(eye(channel.Ms[i], ds[k]), 1)
            end
        end
    elseif aux_params["initial_precoders"] == "white"
        for i = 1:channel.I
            served = served_MS_ids(i, cell_assignment)
            Kc = length(served)

            for k in served
                V[k] = sqrt(Ps[i]/(ds[k]*Kc))*eye(channel.Ms[i], ds[k])
            end
        end
    elseif aux_params["initial_precoders"] == "zeros"
        for i = 1:channel.I
            served = served_MS_ids(i, cell_assignment)
            Kc = length(served)

            for k in served
                V[k] = zeros(channel.Ms[i], ds[k])
            end
        end
    elseif aux_params["initial_precoders"] == "eigendirection"
        for i = 1:channel.I
            served = served_MS_ids(i, cell_assignment)
            Kc = length(served)

            for k in served
                _, _, Vtmp = svd(channel.H[k,i])
                V[k] = sqrt(Ps[i]/Kc)*Vtmp[:, 1:ds[k]]/vecnorm(Vtmp[:, 1:ds[k]])
            end
        end
    end

    return V
end
