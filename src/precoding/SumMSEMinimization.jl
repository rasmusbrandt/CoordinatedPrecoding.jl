immutable SumMSEMinimizationState
    U::Array{Matrix{Complex128},1}
    W::Array{Hermitian{Complex128},1}
    V::Array{Matrix{Complex128},1}
end

function SumMSEMinimization(channel::SinglecarrierChannel, network::Network,
    cell_assignment::CellAssignment, settings=Dict())

    settings = check_and_defaultize_settings(SumMSEMinimizationState, settings)

    Ps = get_transmit_powers(network)
    sigma2s = get_receiver_noise_powers(network)
    ds = get_no_streams(network)

    state = SumMSEMinimizationState(
        Array(Matrix{Complex128}, channel.K),
        unity_MSE_weights(ds),
        initial_precoders(channel, Ps, sigma2s, ds, cell_assignment, settings))
    logdet_rates = Array(Float64, channel.K, maximum(ds), settings["stop_crit"])
    MMSE_rates = Array(Float64, channel.K, maximum(ds), settings["stop_crit"])

    for iter = 1:(settings["stop_crit"]-1)
        update_MSs!(state, channel, sigma2s, cell_assignment)
        logdet_rates[:,:,iter] = calculate_logdet_rates(state)
        MMSE_rates[:,:,iter] = calculate_MMSE_rates(state)
        update_BSs!(state, channel, Ps, cell_assignment, settings)
    end
    update_MSs!(state, channel, sigma2s, cell_assignment)
    logdet_rates[:,:,end] = calculate_logdet_rates(state)
    MMSE_rates[:,:,end] = calculate_MMSE_rates(state)

    if settings["output_protocol"] == 1
        return [ "logdet_rates" => logdet_rates, "MMSE_rates" => MMSE_rates ]
    elseif settings["output_protocol"] == 2
        return [ "logdet_rates" => logdet_rates[:,:,end],
                 "MMSE_rates" => MMSE_rates[:,:,end] ]
    end
end

function check_and_defaultize_settings(::Type{SumMSEMinimizationState}, settings)
    settings = copy(settings)

    # Global settings and consistency checks
    if !haskey(settings, "output_protocol")
        settings["output_protocol"] = 1
    end
    if !haskey(settings, "stop_crit")
        settings["stop_crit"] = 20
    end
    if !haskey(settings, "initial_precoders")
        settings["initial_precoders"] = "dft"
    end
    if settings["output_protocol"] != 1 && settings["output_protocol"] != 2
        error("Unknown output protocol")
    end

    # Local settings
    if !haskey(settings, "SumMSEMinimization:bisection_Gamma_cond")
        settings["SumMSEMinimization:bisection_Gamma_cond"] = 1e10
    end
    if !haskey(settings, "SumMSEMinimization:bisection_singular_Gamma_mu_lower_bound")
        settings["SumMSEMinimization:bisection_singular_Gamma_mu_lower_bound"] = 1e-14
    end
    if !haskey(settings, "SumMSEMinimization:bisection_max_iters")
        settings["SumMSEMinimization:bisection_max_iters"] = 5e1
    end
    if !haskey(settings, "SumMSEMinimization:bisection_tolerance")
        settings["SumMSEMinimization:bisection_tolerance"] = 1e-3
    end

    return settings
end

function update_MSs!(state::SumMSEMinimizationState, channel::SinglecarrierChannel,
    sigma2s::Vector{Float64}, cell_assignment::CellAssignment)

    ds = [ size(state.W[k], 1) for k = 1:channel.K ]

    for i = 1:channel.I
        for k in served_MS_ids(i, cell_assignment)
            # Received covariance
            Phi = Hermitian(complex(sigma2s[k]*eye(channel.Ns[k])))
            for j = 1:channel.I
                for l in served_MS_ids(j, cell_assignment)
                    #Phi += Hermitian(channel.H[k,j]*(state.V[l]*state.V[l]')*channel.H[k,j]')
                    herk!(Phi.uplo, 'N', complex(1.), channel.H[k,j]*state.V[l], complex(1.), Phi.S)
                end
            end

            # MMSE receiver and optimal MSE weight
            F = channel.H[k,i]*state.V[k]
            state.U[k] = Phi\F
            state.W[k] = Hermitian((eye(ds[k]) - state.U[k]'*F)\eye(ds[k]))
        end
    end
end

function update_BSs!(state::SumMSEMinimizationState, channel::SinglecarrierChannel, 
    Ps::Vector{Float64}, cell_assignment::CellAssignment, settings)

    for i = 1:channel.I
        # Virtual uplink covariance
        Gamma = Hermitian(complex(zeros(channel.Ms[i],channel.Ms[i])))
        for j = 1:channel.I
            for l in served_MS_ids(j, cell_assignment)
                Gamma += Hermitian(channel.H[l,i]'*(state.U[l]*state.U[l]')*channel.H[l,i])
            end
        end

        # Find optimal Lagrange multiplier
        mu_star, Gamma_eigen =
            optimal_mu(i, Gamma, state, channel, Ps, cell_assignment, settings)

        # Precoders (reuse EVD)
        for k in served_MS_ids(i, cell_assignment)
            state.V[k] = Gamma_eigen.vectors*Diagonal(1./(Gamma_eigen.values .+ mu_star))*Gamma_eigen.vectors'*channel.H[k,i]'*state.U[k]
        end
    end
end

function optimal_mu(i::Int, Gamma::Hermitian{Complex128},
    state::SumMSEMinimizationState, channel::SinglecarrierChannel,
    Ps::Vector{Float64}, cell_assignment::CellAssignment, settings)

    # Build bisector function
    bis_M = Hermitian(complex(zeros(channel.Ms[i], channel.Ms[i])))
    for k in served_MS_ids(i, cell_assignment)
        #bis_M += Hermitian(channel.H[k,i]'*(state.U[k]*state.U[k]')*channel.H[k,i])
        herk!(bis_M.uplo, 'N', complex(1.), channel.H[k,i]'*state.U[k], complex(1.), bis_M.S)
    end
    Gamma_eigen = eigfact(Gamma)
    bis_JMJ_diag = real(diag(Gamma_eigen.vectors'*bis_M*Gamma_eigen.vectors))
    f(mu) = sum(bis_JMJ_diag./((Gamma_eigen.values .+ mu).*(Gamma_eigen.values .+ mu)))

    # mu lower bound
    if abs(maximum(Gamma_eigen.values))/abs(minimum(Gamma_eigen.values)) < settings["SumMSEMinimization:bisection_Gamma_cond"]
        # Gamma is invertible
        mu_lower = 0
    else
        mu_lower = settings["SumMSEMinimization:bisection_singular_Gamma_mu_lower_bound"]
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
        while no_iters < settings["SumMSEMinimization:bisection_max_iters"]
            conv_crit = (Ps[i] - f(mu_upper))/Ps[i]

            if conv_crit < settings["SumMSEMinimization:bisection_tolerance"]
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

        if no_iters == settings["SumMSEMinimization:bisection_max_iters"]
            println("Power bisection: reached max iterations.")
        end

        # The upper point is always feasible, therefore we use it
        return mu_upper, Gamma_eigen
    end
end
