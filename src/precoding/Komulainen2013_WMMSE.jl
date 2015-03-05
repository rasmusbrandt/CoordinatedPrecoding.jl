immutable Komulainen2013_WMMSEState
    U::Array{Matrix{Complex128},1}
    W::Array{Hermitian{Complex128},1}
    W_diag::Array{Diagonal{Float64},1}
    V::Array{Matrix{Complex128},1}
end

function Komulainen2013_WMMSE(channel, network)
    assignment = get_assignment(network)

    K = get_no_MSs(network)
    Ps = get_transmit_powers(network)
    sigma2s = get_receiver_noise_powers(network)
    ds = get_no_streams(network)

    aux_params = get_aux_precoding_params(network)
    @defaultize_param! aux_params "Komulainen2013_WMMSE:bisection_Gamma_cond" 1e10
    @defaultize_param! aux_params "Komulainen2013_WMMSE:bisection_singular_Gamma_mu_lower_bound" 1e-14
    @defaultize_param! aux_params "Komulainen2013_WMMSE:bisection_max_iters" 5e1
    @defaultize_param! aux_params "Komulainen2013_WMMSE:bisection_tolerance" 1e-3

    state = Komulainen2013_WMMSEState(
        Array(Matrix{Complex128}, K),
        unity_MSE_weights(ds),
        Array(Diagonal{Float64}, K),
        initial_precoders(channel, Ps, sigma2s, ds, assignment, aux_params))
    objective = Float64[]
    logdet_rates = Array(Float64, K, maximum(ds), aux_params["max_iters"])
    MMSE_rates = Array(Float64, K, maximum(ds), aux_params["max_iters"])
    allocated_power = Array(Float64, K, maximum(ds), aux_params["max_iters"])

    iters = 0; conv_crit = Inf
    while iters < aux_params["max_iters"]
        update_MSs!(state, channel, sigma2s, assignment)
        iters += 1

        # Results after this iteration
        logdet_rates[:,:,iters] = calculate_logdet_rates(state)
        push!(objective, sum(logdet_rates[:,:,iters]))
        MMSE_rates[:,:,iters] = calculate_MMSE_rates(state)
        allocated_power[:,:,iters] = calculate_allocated_power(state)

        # Check convergence
        if iters >= 2
            conv_crit = abs(objective[end] - objective[end-1])/abs(objective[end-1])
            if conv_crit < aux_params["stop_crit"]
                Lumberjack.debug("Komulainen2013_WMMSE converged.",
                    [ :no_iters => iters, :final_objective => objective[end],
                      :conv_crit => conv_crit, :stop_crit => aux_params["stop_crit"],
                      :max_iters => aux_params["max_iters"] ])
                break
            end
        end

        # Begin next iteration, unless the loop will end
        if iters < aux_params["max_iters"]
            update_BSs!(state, channel, Ps, assignment, aux_params)
        end
    end
    if iters == aux_params["max_iters"]
        Lumberjack.debug("Komulainen2013_WMMSE did NOT converge.",
            [ :no_iters => iters, :final_objective => objective[end],
              :conv_crit => conv_crit, :stop_crit => aux_params["stop_crit"],
              :max_iters => aux_params["max_iters"] ])
    end

    results = PrecodingResults()
    if aux_params["output_protocol"] == :all_iterations
        results["objective"] = objective
        results["logdet_rates"] = logdet_rates
        results["MMSE_rates"] = MMSE_rates
        results["allocated_power"] = allocated_power
    elseif aux_params["output_protocol"] == :final_iteration
        results["objective"] = objective[iters]
        results["logdet_rates"] = logdet_rates[:,:,iters]
        results["MMSE_rates"] = MMSE_rates[:,:,iters]
        results["allocated_power"] = allocated_power[:,:,iters]
    end
    return results
end

function update_MSs!(state::Komulainen2013_WMMSEState,
    channel::SinglecarrierChannel, sigma2s, assignment)

    ds = [ size(state.W[k], 1) for k = 1:channel.K ]

    for i = 1:channel.I; for k in served_MS_ids(i, assignment)
        Phi = Hermitian(complex(sigma2s[k]*eye(channel.Ns[k])))
        for j = 1:channel.I; for l in served_MS_ids(j, assignment)
            #Phi += Hermitian(channel.H[k,j]*(state.V[l]*state.V[l]')*channel.H[k,j]')
            Base.LinAlg.BLAS.herk!(Phi.uplo, 'N', complex(1.), channel.H[k,j]*state.V[l], complex(1.), Phi.S)
        end; end

        # MMSE receiver, optimal MSE weight, and actually used (diagonal) MSE weight
        F = channel.H[k,i]*state.V[k]
        state.U[k] = Phi\F
        Emmse = eye(ds[k]) - state.U[k]'*F
        state.W[k] = Hermitian(Emmse\eye(ds[k]))
        state.W_diag[k] = Diagonal(max(1, abs(1./diag(Emmse))))
    end; end
end

function update_BSs!(state::Komulainen2013_WMMSEState,
    channel::SinglecarrierChannel, Ps, assignment, aux_params)

    for i in active_BSs(assignment)
        Gamma = Hermitian(zeros(Complex128, channel.Ms[i], channel.Ms[i]))
        for j = 1:channel.I; for l in served_MS_ids(j, assignment)
            Gamma += Hermitian(channel.H[l,i]'*(state.U[l]*state.W_diag[l]*state.U[l]')*channel.H[l,i])
        end; end

        # Find optimal Lagrange multiplier
        mu_star, Gamma_eigen =
            optimal_mu(i, Gamma, state, channel, Ps, assignment, aux_params)

        # Precoders (reuse EVD)
        for k in served_MS_ids(i, assignment)
            state.V[k] = Gamma_eigen.vectors*Diagonal(1./(abs(Gamma_eigen.values) .+ mu_star))*Gamma_eigen.vectors'*channel.H[k,i]'*state.U[k]*state.W_diag[k]
        end
    end
end

function optimal_mu(i, Gamma, state::Komulainen2013_WMMSEState,
    channel::SinglecarrierChannel, Ps, assignment, aux_params)

    # Build bisector function
    bis_M = Hermitian(zeros(Complex128, channel.Ms[i], channel.Ms[i]))
    for k in served_MS_ids(i, assignment)
        #bis_M += Hermitian(channel.H[k,i]'*(state.U[k]*(state.W_diag[k]*state.W_diag[k])*state.U[k]')*channel.H[k,i])
        Base.LinAlg.BLAS.herk!(bis_M.uplo, 'N', complex(1.), channel.H[k,i]'*state.U[k]*state.W_diag[k], complex(1.), bis_M.S)
    end
    Gamma_eigen = eigfact(Gamma); Gamma_eigen_values = abs(Gamma_eigen.values)
    bis_JMJ_diag = abs(diag(Gamma_eigen.vectors'*bis_M*Gamma_eigen.vectors))
    f(mu) = sum(bis_JMJ_diag./((Gamma_eigen_values .+ mu).*(Gamma_eigen_values .+ mu)))

    # mu lower bound
    if maximum(Gamma_eigen_values)/minimum(Gamma_eigen_values) < aux_params["Komulainen2013_WMMSE:bisection_Gamma_cond"]
        # Gamma is invertible
        mu_lower = 0.
    else
        mu_lower = aux_params["Komulainen2013_WMMSE:bisection_singular_Gamma_mu_lower_bound"]
    end

    if f(mu_lower) <= Ps[i]
        # No bisection needed
        return mu_lower, Gamma_eigen
    else
        # mu upper bound
        mu_upper = sqrt(channel.Ms[i]/Ps[i]*maximum(bis_JMJ_diag)) - minimum(Gamma_eigen_values)
        if f(mu_upper) > Ps[i]
            Lumberjack.error("Power bisection: infeasible mu upper bound.")
        end

        no_iters = 0
        while no_iters < aux_params["Komulainen2013_WMMSE:bisection_max_iters"]
            conv_crit = (Ps[i] - f(mu_upper))/Ps[i]

            if conv_crit < aux_params["Komulainen2013_WMMSE:bisection_tolerance"]
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

        if no_iters == aux_params["Komulainen2013_WMMSE:bisection_max_iters"]
            Lumberjack.warn("Power bisection: reached max iterations.")
        end

        # The upper point is always feasible, therefore we use it
        return mu_upper, Gamma_eigen
    end
end
