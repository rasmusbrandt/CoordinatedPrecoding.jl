immutable Gomadam2008_MinWLIState
    U::Array{Matrix{Complex128},1}
    W::Array{Hermitian{Complex128},1} # these are only used for rate calculations
    V::Array{Matrix{Complex128},1}
end

function Gomadam2008_MinWLI(channel, network)
    assignment = get_assignment(network)

    K = get_no_MSs(network)
    Ps = get_transmit_powers(network)
    sigma2s = get_receiver_noise_powers(network)
    ds = get_no_streams(network); max_d = maximum(ds)
    alphas = get_user_priorities(network)
    aux_params = get_aux_precoding_params(network)

    state = Gomadam2008_MinWLIState(
        Array(Matrix{Complex128}, K),
        Array(Hermitian{Complex128}, K),
        initial_precoders(channel, Ps, sigma2s, ds, assignment, aux_params))
    objective = Float64[]
    logdet_rates = Array(Float64, K, max_d, aux_params["max_iters"])
    MMSE_rates = Array(Float64, K, max_d, aux_params["max_iters"])
    weighted_logdet_rates = Array(Float64, K, max_d, aux_params["max_iters"])
    weighted_MMSE_rates = Array(Float64, K, max_d, aux_params["max_iters"])
    allocated_power = Array(Float64, K, max_d, aux_params["max_iters"])

    iters = 0; conv_crit = Inf
    while iters < aux_params["max_iters"]
        update_MSs!(state, channel, sigma2s, assignment)
        iters += 1

        # Results after this iteration
        logdet_rates[:,:,iters] = calculate_logdet_rates(state)
        push!(objective, sum(logdet_rates[:,:,iters]))
        MMSE_rates[:,:,iters] = calculate_MMSE_rates(state)
        weighted_logdet_rates[:,:,iters] = calculate_weighted_logdet_rates(state, alphas)
        weighted_MMSE_rates[:,:,iters] = calculate_weighted_MMSE_rates(state, alphas)
        allocated_power[:,:,iters] = calculate_allocated_power(state)

        # Check convergence
        if iters >= 2
            conv_crit = abs(objective[end] - objective[end-1])/abs(objective[end-1])
            if conv_crit < aux_params["stop_crit"]
                Lumberjack.debug("Gomadam2008_MinWLI converged.",
                    [ :no_iters => iters,
                      :final_objective => objective[end],
                      :conv_crit => conv_crit,
                      :stop_crit => aux_params["stop_crit"],
                      :max_iters => aux_params["max_iters"] ]
                )
                break
            end
        end

        # Begin next iteration, unless the loop will end
        if iters < aux_params["max_iters"]
            update_BSs!(state, channel, Ps, sigma2s, assignment, aux_params)
        end
    end
    if iters == aux_params["max_iters"]
        Lumberjack.debug("Gomadam2008_MinWLI did NOT converge.",
            [ :no_iters => iters,
              :final_objective => objective[end],
              :conv_crit => conv_crit,
              :stop_crit => aux_params["stop_crit"],
              :max_iters => aux_params["max_iters"] ]
        )
    end

    results = PrecodingResults()
    if aux_params["output_protocol"] == :all_iterations
        results["objective"] = objective
        results["logdet_rates"] = logdet_rates
        results["MMSE_rates"] = MMSE_rates
        results["weighted_logdet_rates"] = weighted_logdet_rates
        results["weighted_MMSE_rates"] = weighted_MMSE_rates
        results["allocated_power"] = allocated_power
    elseif aux_params["output_protocol"] == :final_iteration
        results["objective"] = objective[iters]
        results["logdet_rates"] = logdet_rates[:,:,iters]
        results["MMSE_rates"] = MMSE_rates[:,:,iters]
        results["weighted_logdet_rates"] = weighted_logdet_rates[:,:,iters]
        results["weighted_MMSE_rates"] = weighted_MMSE_rates[:,:,iters]
        results["allocated_power"] = allocated_power[:,:,iters]
    end
    return results
end

function update_MSs!(state::Gomadam2008_MinWLIState,
    channel::SinglecarrierChannel, sigma2s, assignment)

    ds = [ size(state.V[k], 2) for k = 1:channel.K ]

    for i = 1:channel.I; for k in served_MS_ids(i, assignment)
        Psi = zeros(Complex128, channel.Ns[k], channel.Ns[k])
        for j = 1:channel.I; for l in served_MS_ids_except_me(k, j, assignment)
            F = channel.H[k,j]*state.V[l]
            Psi += F*F'
        end; end

        # Zero-forcing receiver
        state.U[k] = eigfact(Hermitian(Psi), 1:ds[k]).vectors

        # Optimal MSE weights
        F = channel.H[k,i]*state.V[k]
        Ummse = Hermitian(Psi + F*F' + sigma2s[k]*eye(Psi))\F
        state.W[k] = inv(Hermitian(eye(ds[k]) - Ummse'*F))
    end; end
end

function update_BSs!(state::Gomadam2008_MinWLIState,
    channel::SinglecarrierChannel, Ps, sigma2s, assignment, aux_params)

    ds = [ size(state.W[k], 1) for k = 1:channel.K ]

    for i in active_BSs(assignment)
        Gamma = zeros(Complex128, channel.Ms[i], channel.Ms[i])
        for j = 1:channel.I; for l = served_MS_ids(j, assignment)
            G = channel.H[l,i]'*state.U[l]
            Gamma += G*G'
        end; end

        # Precoders for all served users
        served = served_MS_ids(i, assignment)
        Nserved = length(served)
        for k in served
            G = channel.H[k,i]'*state.U[k]

            # Precoder
            state.V[k] = sqrt(Ps[i]/(Nserved*ds[k]))*eigfact(Hermitian(Gamma - G*G'), 1:ds[k]).vectors
        end
    end
end
