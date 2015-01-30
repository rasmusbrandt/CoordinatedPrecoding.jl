immutable Gomadam2008_MaxSINRState
    U::Array{Matrix{Complex128},1}
    W::Array{Hermitian{Complex128},1} # these are only used for rate calculations
    V::Array{Matrix{Complex128},1}
end

function Gomadam2008_MaxSINR(channel::SinglecarrierChannel, network::Network,
    cell_assignment::CellAssignment, settings=Dict())

    check_and_defaultize_precoding_settings!(settings, Gomadam2008_MaxSINRState)

    K = get_no_MSs(network)
    Ps = get_transmit_powers(network)
    sigma2s = get_receiver_noise_powers(network)
    ds = get_no_streams(network)

    state = Gomadam2008_MaxSINRState(
        zero_receivers(channel, ds),
        unity_MSE_weights(ds),
        initial_precoders(channel, Ps, sigma2s, ds, cell_assignment, settings))
    objective = Float64[]
    logdet_rates = Array(Float64, K, maximum(ds), settings["max_iters"])
    MMSE_rates = Array(Float64, K, maximum(ds), settings["max_iters"])
    allocated_power = Array(Float64, K, maximum(ds), settings["max_iters"])

    iters = 0; conv_crit = Inf
    while iters < settings["max_iters"]
        update_MSs!(state, channel, sigma2s, cell_assignment)
        iters += 1

        # Results after this iteration
        logdet_rates[:,:,iters], t = calculate_logdet_rates(state, settings)
        push!(objective, t)
        MMSE_rates[:,:,iters], _ = calculate_MMSE_rates(state, settings)
        allocated_power[:,:,iters] = calculate_allocated_power(state)

        # Check convergence
        if iters >= 2
            conv_crit = abs(objective[end] - objective[end-1])/abs(objective[end-1])
            if conv_crit < settings["stop_crit"]
                Lumberjack.debug("Gomadam2008_MaxSINR converged.",
                    { :no_iters => iters, :final_objective => objective[end],
                      :conv_crit => conv_crit, :stop_crit => settings["stop_crit"],
                      :max_iters => settings["max_iters"] })
                break
            end
        end

        # Begin next iteration, unless the loop will end
        if iters < settings["max_iters"]
            update_BSs!(state, channel, Ps, sigma2s, cell_assignment, settings)
        end
    end
    if iters == settings["max_iters"]
        Lumberjack.debug("Gomadam2008_MaxSINR did NOT converge.",
            { :no_iters => iters, :final_objective => objective[end],
              :conv_crit => conv_crit, :stop_crit => settings["stop_crit"],
              :max_iters => settings["max_iters"] })
    end

    results = Dict{ASCIIString, Any}()
    if settings["output_protocol"] == 1
        results["objective"] = objective
        results["logdet_rates"] = logdet_rates
        results["MMSE_rates"] = MMSE_rates
        results["allocated_power"] = allocated_power
    elseif settings["output_protocol"] == 2
        results["objective"] = objective[iters]
        results["logdet_rates"] = logdet_rates[:,:,iters]
        results["MMSE_rates"] = MMSE_rates[:,:,iters]
        results["allocated_power"] = allocated_power[:,:,iters]
    end
    return results
end

function check_and_defaultize_precoding_settings!(settings, ::Type{Gomadam2008_MaxSINRState})
    # Global settings
    check_and_defaultize_precoding_settings!(settings)

    # Local settings
end

function update_MSs!(state::Gomadam2008_MaxSINRState,
    channel::SinglecarrierChannel, sigma2s::Vector{Float64},
    cell_assignment::CellAssignment)

    ds = [ size(state.W[k], 1) for k = 1:channel.K ]

    for i = 1:channel.I
        for k in served_MS_ids(i, cell_assignment)
            Phi = Hermitian(complex(sigma2s[k]*eye(channel.Ns[k])))
            for j = 1:channel.I
                for l in served_MS_ids(j, cell_assignment)
                    #Phi += Hermitian(channel.H[k,j]*(state.V[l]*state.V[l]')*channel.H[k,j]')
                    Base.LinAlg.BLAS.herk!(Phi.uplo, 'N', complex(1.), channel.H[k,j]*state.V[l], complex(1.), Phi.S)
                end
            end

            # Per-stream receivers
            for n = 1:ds[k]
                Phi_i_plus_n = Hermitian(
                    Base.LinAlg.BLAS.herk!(Phi.uplo, 'N', complex(-1.), channel.H[k,i]*state.V[k][:,n], complex(1.), copy(Phi.S)),
                    Phi.uplo)
                u = Phi_i_plus_n\channel.H[k,i]*state.V[k][:,n]
                state.U[k][:,n] = u/norm(u,2)
            end

            # Optimal MSE weights (for rate calculation only)
            F = channel.H[k,i]*state.V[k]
            Ummse = Phi\F
            state.W[k] = Hermitian((eye(ds[k]) - Ummse'*F)\eye(ds[k]))
        end
    end
end

function update_BSs!(state::Gomadam2008_MaxSINRState,
    channel::SinglecarrierChannel, Ps::Vector{Float64},
    sigma2s::Vector{Float64}, cell_assignment::CellAssignment, settings)

    ds = [ size(state.W[k], 1) for k = 1:channel.K ]

    for i = 1:channel.I
        # Virtual uplink covariance
        Gamma = Hermitian(complex(zeros(channel.Ms[i],channel.Ms[i])))
        for j = 1:channel.I
            for l = served_MS_ids(j, cell_assignment)
                #Gamma += Hermitian(channel.H[k,i]'*(state.U[k]*state.U[k]')*channel.H[k,i])
                Base.LinAlg.BLAS.herk!(Gamma.uplo, 'N', complex(1.), channel.H[l,i]'*state.U[l], complex(1.), Gamma.S)
            end
        end

        # Per-stream precoders
        served = served_MS_ids(i, cell_assignment)
        Nserved = length(served)
        for k in served
            for n = 1:ds[k]
                Gamma_i_plus_n = Hermitian(
                    Base.LinAlg.BLAS.herk!(Gamma.uplo, 'N', complex(-1.), channel.H[k,i]'*state.U[k][:,n], complex(1.), copy(Gamma.S)) + (sigma2s[k]/Ps[i])*eye(channel.Ms[i]),
                    Gamma.uplo)
                v = Gamma_i_plus_n\channel.H[k,i]'*state.U[k][:,n]
                state.V[k][:,n] = sqrt(Ps[i]/(Nserved*ds[k]))*v/norm(v,2)
            end
        end
    end
end
