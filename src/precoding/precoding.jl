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
include("Gomadam2008_MinWLI.jl")
include("Komulainen2013_WMMSE.jl")
include("Razaviyayn2013_MinMaxWMMSE.jl")
include("Shi2011_WMMSE.jl")
include("SumMSEMinimization.jl")

##########################################################################
# Standard functions to calculate rates from optimal MSE weights
function calculate_logdet_rates(state)
    K = length(state.W)
    ds = Int[ size(state.W[k], 1) for k = 1:K ]; max_d = maximum(ds)

    logdet_rates = zeros(Float64, K, max_d)
    for k = 1:K; if ds[k] > 0
        # W is p.d., so we should only get abs eigenvalues. Numerically we may
        # get some imaginary noise however. Also, numerically the eigenvalues
        # may be less than 1, so we need to handle that to not get negative
        # rates.
        r = log2(max(1, abs(eigvals(state.W[k]))))

        if ds[k] < max_d
            logdet_rates[k,:] = cat(1, r, zeros(Float64, max_d - ds[k]))
        else
            logdet_rates[k,:] = r
        end
    end; end

    return logdet_rates
end

function calculate_MMSE_rates(state)
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
        r = log2(max(1, abs(1./diag(E))))

        if ds[k] < max_d
            MMSE_rates[k,:] = cat(1, r, zeros(Float64, max_d - ds[k]))
        else
            MMSE_rates[k,:] = r
        end
    end; end

    return MMSE_rates
end

function calculate_allocated_power(state)
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
function initial_receivers(channel::SinglecarrierChannel, Ps, sigma2s, ds,
    assignment, aux_params)

    U = Array(Matrix{Complex128}, channel.K)

    if haskey(aux_params, "initial_receivers")
        if aux_params["initial_receivers"] == "eigendirection"
            for i = 1:channel.I; for k in served_MS_ids(i, assignment)
                Utmp, _, _ = svd(channel.H[k,i])
                U[k] = Utmp[:, 1:ds[k]]
            end; end

            return U
        elseif aux_params["initial_receivers"] == "white"
            for i = 1:channel.I; for k in served_MS_ids(i, assignment)
                U[k] = eye(channel.Ns[k], ds[k])
            end; end

            return U
        elseif aux_params["initial_receivers"] == "zeros"
            for i = 1:channel.I; for k in served_MS_ids(i, assignment)
                U[k] = zeros(channel.Ns[k], ds[k])
            end; end

            return U
        end
    end

    # Default: dft precoders
    for i = 1:channel.I; for k in served_MS_ids(i, assignment)
        U[k] = fft(eye(channel.Ns[k], ds[k]), 1)
    end; end

    return U
end

function initial_MSE_weights(channel::SinglecarrierChannel, Ps, sigma2s, ds,
    assignment, aux_params)

    # So far only the trivial initial point available
    return [ Hermitian(eye(ds[k])) for k = 1:channel.K ]
end

function initial_precoders(channel::SinglecarrierChannel, Ps, sigma2s, ds,
    assignment, aux_params)

    V = Array(Matrix{Complex128}, channel.K)

    if haskey(aux_params, "initial_precoders")
        if aux_params["initial_precoders"] == "eigendirection"
            for i = 1:channel.I
                served = served_MS_ids(i, assignment)
                Kc = length(served)

                for k in served
                    _, _, Vtmp = svd(channel.H[k,i])
                    V[k] = sqrt(Ps[i]/Kc)*Vtmp[:, 1:ds[k]]/vecnorm(Vtmp[:, 1:ds[k]])
                end
            end

            return V
        elseif aux_params["initial_precoders"] == "white"
            for i = 1:channel.I
                served = served_MS_ids(i, assignment)
                Kc = length(served)

                for k in served
                    V[k] = sqrt(Ps[i]/(ds[k]*Kc))*eye(channel.Ms[i], ds[k])
                end
            end

            return V
        elseif aux_params["initial_precoders"] == "zeros"
            for i = 1:channel.I
                served = served_MS_ids(i, assignment)
                Kc = length(served)

                for k in served
                    V[k] = zeros(channel.Ms[i], ds[k])
                end
            end

            return V
        end
    end

    # Default: dft precoders
    for i = 1:channel.I
        served = served_MS_ids(i, assignment)
        Kc = length(served)

        for k in served
            V[k] = sqrt(Ps[i]/(channel.Ms[i]*ds[k]*Kc))*fft(eye(channel.Ms[i], ds[k]), 1)
        end
    end

    return V
end
