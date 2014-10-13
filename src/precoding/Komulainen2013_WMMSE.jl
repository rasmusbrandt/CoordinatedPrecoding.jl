immutable Komulainen2013_WMMSEState
    A::Array{Matrix{Complex128},1}
    W::Matrix{Float64}
    V::Array{Matrix{Complex128},1}

    Phi::Array{Hermitian{Complex128},1}
    Gamma::Array{Hermitian{Complex128},1}
end

function Komulainen2013_WMMSE(channel::SinglecarrierChannel, network::Network,
    cell_assignment::CellAssignment, settings=Dict())

    Ps = get_transmit_powers(network)
    sigma2s = get_receiver_noise_powers(network)
    ds = get_no_streams(network)

    defaultize_settings!(Komulainen2013_WMMSEState, settings)

    state = Komulainen2013_WMMSEState(
        Array(Matrix{Complex128}, channel.K),
        ones(Float64, channel.K, maximum(ds)), # Unused MSE weights (i.e. weights that are never updated) are kept at unity.
        initial_precoders(channel, Ps, sigma2s, ds, cell_assignment, settings),
        Array(Hermitian{Complex128}, channel.K),
        Array(Hermitian{Complex128}, channel.K))
    user_rates = Array(Float64, channel.K, maximum(ds), settings["stop_crit"])

    for iter = 1:(settings["stop_crit"]-1)
        update_MSs!(state, channel, sigma2s, ds, cell_assignment)
        user_rates[:,:,iter] = calculate_user_rates(state)
        update_BSs!(state, channel, Ps, cell_assignment, settings)
    end
    update_MSs!(state, channel, sigma2s, ds, cell_assignment)
    user_rates[:,:,end] = calculate_user_rates(state)

    return user_rates
end

function defaultize_settings!(::Type{Komulainen2013_WMMSEState}, settings)
    if !haskey(settings, "stop_crit")
        settings["stop_crit"] = 20
    end
    if !haskey(settings, "initial_precoders")
        settings["initial_precoders"] = "dft"
    end
    if !haskey(settings, "bisection_Gamma_cond")
        settings["bisection_Gamma_cond"] = 1e10
    end
    if !haskey(settings, "bisection_singular_Gamma_mu_lower_bound")
        settings["bisection_singular_Gamma_mu_lower_bound"] = 1e-14
    end
    if !haskey(settings, "bisection_max_iters")
        settings["bisection_max_iters"] = 5e1
    end
    if !haskey(settings, "bisection_tolerance")
        settings["bisection_tolerance"] = 1e-3
    end
end

function update_MSs!(state::Komulainen2013_WMMSEState, channel::SinglecarrierChannel,
    sigma2s::Vector{Float64}, ds::Vector{Int},
    cell_assignment::CellAssignment)

    for k = 1:channel.K
        # Received covariance
        state.Phi[k] = Hermitian(complex(sigma2s[k]*eye(channel.Ns[k])))
        for j = 1:channel.I
            for l in served_MS_ids(j, cell_assignment)
                #state.Phi[k] += Hermitian(channel.H[k,j]*(state.V[l]*state.V[l]')*channel.H[k,j]')
                herk!(state.Phi[k].uplo, 'N', complex(1.), channel.H[k,j]*state.V[l], complex(1.), state.Phi[k].S)
            end
        end

        # MMSE receiver and optimal MSE weight
        i = serving_BS_id(k, cell_assignment)
        Fk = channel.H[k,i]*state.V[k]
        state.A[k] = state.Phi[k]\Fk
        for n = 1:ds[k]
            state.W[k,n] = abs(1/(1 - state.A[k][:,n]'*Fk[:,n])[1])
        end
    end
end

function update_BSs!(state::Komulainen2013_WMMSEState, channel::SinglecarrierChannel, 
    Ps::Vector{Float64}, cell_assignment::CellAssignment, settings)

    for i = 1:channel.I
        # Covariance
        state.Gamma[i] = Hermitian(complex(zeros(channel.Ms[i],channel.Ms[i])))
        for k = 1:channel.K
            Wk = diagm(squeeze(state.W[k,1:size(state.A[k],2)], 1))
            state.Gamma[i] += Hermitian(channel.H[k,i]'*(state.A[k]*Wk*state.A[k]')*channel.H[k,i])
        end

        # Find optimal Lagrange multiplier
        mu_star, Gamma_eigen =
            optimal_mu(i, state, channel, Ps, cell_assignment, settings)

        # Precoders (reuse EVD)
        for k in served_MS_ids(i, cell_assignment)
            Wk = diagm(squeeze(state.W[k,1:size(state.A[k],2)], 1))
            state.V[k] = Gamma_eigen.vectors*Diagonal(1./(Gamma_eigen.values .+ mu_star))*Gamma_eigen.vectors'*channel.H[k,i]'*state.A[k]*Wk
        end
    end
end

function optimal_mu(i::Int, state::Komulainen2013_WMMSEState,
    channel::SinglecarrierChannel, Ps::Vector{Float64},
    cell_assignment::CellAssignment, settings)

    # Build bisector function
    bis_M = Hermitian(complex(zeros(channel.Ms[i], channel.Ms[i])))
    for k in served_MS_ids(i, cell_assignment)
        Wk = diagm(squeeze(state.W[k,1:size(state.A[k],2)], 1))

        #bis_M += Hermitian(channel.H[k,i]'*(state.A[k]*(Wk*Wk)*state.A[k]')*channel.H[k,i])
        herk!(bis_M.uplo, 'N', complex(1.), channel.H[k,i]'*state.A[k]*Wk, complex(1.), bis_M.S)
    end
    Gamma_eigen = eigfact(state.Gamma[i])
    bis_JMJ_diag = real(diag(Gamma_eigen.vectors'*bis_M*Gamma_eigen.vectors))
    f(mu) = sum(bis_JMJ_diag./((Gamma_eigen.values .+ mu).*(Gamma_eigen.values .+ mu)))

    # mu lower bound
    if abs(maximum(Gamma_eigen.values))/abs(minimum(Gamma_eigen.values)) < settings["bisection_Gamma_cond"]
        # Gamma is invertible
        mu_lower = 0
    else
        mu_lower = settings["bisection_singular_Gamma_mu_lower_bound"]
    end

    if f(mu_lower) <= Ps[i]
        # No bisection needed
        return mu_lower, Gamma_eigen
    else
        # mu upper bound
        mu_upper = sqrt(channel.Ms[i]/Ps[i]*maximum(bis_JMJ_diag)) - abs(minimum(Gamma_eigen.values))
        if f(mu_upper) > Ps[i]
            error("Power bisection: infeasible mu upper bound.")
        end

        no_iters = 0
        while no_iters < settings["bisection_max_iters"]
            conv_crit = (Ps[i] - f(mu_upper))/Ps[i]

            if conv_crit < settings["bisection_tolerance"]
                break
            else
                mu = (1/2)*(mu_lower + mu_upper)

                if f(mu) < Ps[i]
                    # New point feasible, replace upper point
                    mu_upper = mu
                else
                    # New point not feasible, replace lower point
                    mu_lower = mu
                end
            end

            no_iters += 1
        end

        if no_iters == settings["bisection_max_iters"]
            println("Power bisection: reached max iterations.")
        end

        # The upper point is always feasible, therefore we use it
        return mu_upper, Gamma_eigen
    end
end

function calculate_user_rates(state::Komulainen2013_WMMSEState)
    K, max_d = size(state.W)

    user_rates = Array(Float64, K, max_d)

    for k = 1:K
        for n = 1:max_d
            user_rates[k,n] = log2(abs(state.W[k,n]))
        end
    end

    return user_rates
end
