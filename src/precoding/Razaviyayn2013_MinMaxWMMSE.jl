if Pkg.installed("Gurobi") != nothing
    import Gurobi
end

immutable Razaviyayn2013_MinMaxWMMSEState
    U::Array{Matrix{Complex128},1}
    W::Array{Hermitian{Complex128},1}
    V::Array{Matrix{Complex128},1}
end

function Razaviyayn2013_MinMaxWMMSE(channel, network)
    assignment = get_assignment(network)

    # The implementation is currently limited in the respects below. This is in
    # order to simplify the Gurobi optimization variable indexing. With equal
    # number of antennas and streams, it is very easy to calculate the variable
    # offsets required. (See functions v_ind and mu_ind.)
    require_equal_num_BS_antennas(network)
    require_equal_num_streams(network)
    require_equal_num_MSs_per_cell(assignment)

    K = get_num_MSs(network)
    Ps = get_transmit_powers(network)
    sigma2s = get_receiver_noise_powers(network)
    ds = get_num_streams(network); max_d = maximum(ds)
    alphas = get_user_priorities(network)

    aux_params = get_aux_precoding_params(network)
    @defaultize_param! aux_params "Razaviyayn2013_MinMaxWMMSE:Gurobi_verbose" false
    @defaultize_param! aux_params "Razaviyayn2013_MinMaxWMMSE:Gurobi_PSDTol" 1e-5

    state = Razaviyayn2013_MinMaxWMMSEState(
        Array(Matrix{Complex128}, K),
        Array(Hermitian{Complex128}, K),
        initial_precoders(channel, Ps, sigma2s, ds, assignment, aux_params))
    objective = Float64[]
    logdet_rates = zeros(Float64, K, max_d, aux_params["max_iters"])
    MMSE_rates = zeros(Float64, K, max_d, aux_params["max_iters"])
    weighted_logdet_rates = zeros(Float64, K, max_d, aux_params["max_iters"])
    weighted_MMSE_rates = zeros(Float64, K, max_d, aux_params["max_iters"])
    allocated_power = zeros(Float64, K, max_d, aux_params["max_iters"])

    if Pkg.installed("Gurobi") == nothing
        iters = 1; objective = 0.
        Lumberjack.warn("Razaviyayn2013_MinMaxWMMSE exiting because Gurobi.jl is not installed.")
    else
        iters = 0; conv_crit = Inf
        while iters < aux_params["max_iters"]
            update_MSs!(state, channel, sigma2s, assignment)
            iters += 1

            # Results after this iteration
            logdet_rates[:,:,iters] = calculate_logdet_rates(state)
            push!(objective, minimum(sum(logdet_rates[:,:,iters], 2)))
            MMSE_rates[:,:,iters] = calculate_MMSE_rates(state)
            weighted_logdet_rates[:,:,iters] = calculate_weighted_logdet_rates(state, alphas)
            weighted_MMSE_rates[:,:,iters] = calculate_weighted_MMSE_rates(state, alphas)
            allocated_power[:,:,iters] = calculate_allocated_power(state)

            # Check convergence
            if iters >= 2
                conv_crit = abs(objective[end] - objective[end-1])/abs(objective[end-1])
                if conv_crit < aux_params["stop_crit"]
                    Lumberjack.debug("Razaviyayn2013_MinMaxWMMSE converged.",
                        @compat Dict(
                            :num_iters => iters,
                            :final_objective => objective[end],
                            :conv_crit => conv_crit,
                            :stop_crit => aux_params["stop_crit"],
                            :max_iters => aux_params["max_iters"]))
                    break
                end
            end

            # Begin next iteration, unless the loop will end
            if iters < aux_params["max_iters"]
                update_BSs!(state, channel, Ps, sigma2s, assignment, aux_params)
            end
        end
        if iters == aux_params["max_iters"]
            Lumberjack.debug("Razaviyayn2013_MinMaxWMMSE did NOT converge.",
                @compat Dict(
                    :num_iters => iters,
                    :final_objective => objective[end],
                    :conv_crit => conv_crit,
                    :stop_crit => aux_params["stop_crit"],
                    :max_iters => aux_params["max_iters"]))
        end
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

function update_MSs!(state::Razaviyayn2013_MinMaxWMMSEState,
    channel::SinglecarrierChannel, sigma2s, assignment)

    ds = [ size(state.V[k], 2) for k = 1:channel.K ]

    for i = 1:channel.I; for k in served_MS_ids(i, assignment)
        Phi = complex(sigma2s[k]*eye(channel.Ns[k]))
        for j = 1:channel.I; for l in served_MS_ids(j, assignment)
            F = channel.H[k,j]*state.V[l]
            Phi += F*F'
        end; end

        # MMSE receiver and optimal MSE weight
        F = channel.H[k,i]*state.V[k]
        state.U[k] = Hermitian(Phi)\F
        Emmse = UniformScaling(1.) - state.U[k]'*F
        state.W[k] = inv(Hermitian((Emmse + Emmse')/2))
    end; end
end

# Some of the vectors describing the power constraints and rate constraints
# could be prebuilt and reused in the iterations. Probably no major performance
# boost would come from that though, based on experience with the
# JointPrecodingMCSSelection implementation.
function update_BSs!(state::Razaviyayn2013_MinMaxWMMSEState,
    channel::SinglecarrierChannel, Ps, sigma2s, assignment, aux_params)

    # The following is OK, since we have checked this at the top of
    # the outer function.
    K = channel.K; I = channel.I; Kc = @compat Int(channel.K/channel.I)
    N = channel.Ns[1]; M = channel.Ms[1]; d = size(state.W[1], 2)

    # Bookkeeping for the Gurobi model
    num_t_vars = 1
    num_v_vars = 2*K*M*d
    num_lin_vars_per_rate_constraint = 2*M*d + 1 # including t variable
    num_quad_vars_per_rate_constraint = K*d*(2*M)^2
    num_vars_per_power_constraint = 2*Kc*M*d

    # Index functions for optimization variables
    t_ind = 1
    v_ind(k, m, n, imag::Bool) = (num_t_vars + (k-1)*M*d*2 + (m-1)*d*2 + (n-1)*2 + 1 + @compat Int(imag))

    # Gurobi environment
    env = Gurobi.Env()
    if !aux_params["Razaviyayn2013_MinMaxWMMSE:Gurobi_verbose"]
        Gurobi.setparam!(env, "OutputFlag", 0)
    end
    Gurobi.setparam!(env, "PSDTol", aux_params["Razaviyayn2013_MinMaxWMMSE:Gurobi_PSDTol"])
    model = Gurobi.Model(env, "Razaviyayn2013_MinMaxWMMSE", finalize_env=true)
    Gurobi.set_sense!(model, :maximize)

    # Objective
    Gurobi.add_cvar!(model, 1., 0, Inf) # t

    # Other variables not in objective
    Gurobi.add_cvars!(model, zeros(Float64, num_v_vars))

    # Need to update model before adding constraints
    Gurobi.update_model!(model)

    # Rate constraints
    rate_constraint_lind = Array(Int, num_lin_vars_per_rate_constraint)
    rate_constraint_lval = Array(Float64, num_lin_vars_per_rate_constraint)
    rate_constraint_qr = Array(Int, num_quad_vars_per_rate_constraint)
    rate_constraint_qc = Array(Int, num_quad_vars_per_rate_constraint)
    rate_constraint_qv = Array(Float64, num_quad_vars_per_rate_constraint)
    for i = 1:I; for k = served_MS_ids(i, assignment)
        # Linear part
        g = (channel.H[k,i]'*state.U[k])*state.W[k]
        l_ind = 1
        for m = 1:M; for n = 1:d
            real_ind = v_ind(k, m, n, false)
            rate_constraint_lind[l_ind] = real_ind
            rate_constraint_lval[l_ind] = -2.0*real(g[m,n])
            l_ind += 1

            imag_ind = v_ind(k, m, n, true)
            rate_constraint_lind[l_ind] = imag_ind
            rate_constraint_lval[l_ind] = -2.0*imag(g[m,n])
            l_ind += 1
        end; end
        rate_constraint_lind[l_ind] = t_ind
        rate_constraint_lval[l_ind] = 1.

        # Quadratic part
        # Enforce PSDness explicitly with get_Gr and get_Gi.
        # Hopefully these are inlined.
        q_ind = 1
        for j = 1:I
            G = (channel.H[k,j]'*state.U[k])*state.W[k]*(state.U[k]'*channel.H[k,j])

            # symmetric
            Gr = real(G)
            get_Gr(mm1, mm2) = ( mm1 >= mm2 ? Gr[mm1,mm2] :  Gr[mm2,mm1] )

            # skew-symmetric
            Gi = imag(G)
            get_Gi(mm1, mm2) = ( mm1 >= mm2 ? Gi[mm1,mm2] : -Gi[mm2,mm1] )

            for l = served_MS_ids(j, assignment)
                for nn = 1:d
                    for mm1 = 1:M
                        real_ind1 = v_ind(l, mm1, nn, false)
                        imag_ind1 = v_ind(l, mm1, nn, true)
                        for mm2 = 1:M
                            real_ind2 = v_ind(l, mm2, nn, false)
                            imag_ind2 = v_ind(l, mm2, nn, true)

                            # real(v)'*real(Q)*real(v)
                            rate_constraint_qr[q_ind] = real_ind1
                            rate_constraint_qc[q_ind] = real_ind2
                            rate_constraint_qv[q_ind] = get_Gr(mm1, mm2)
                            q_ind += 1

                            # imag(v)'*real(Q)*imag(v)
                            rate_constraint_qr[q_ind] = imag_ind1
                            rate_constraint_qc[q_ind] = imag_ind2
                            rate_constraint_qv[q_ind] = get_Gr(mm1, mm2)
                            q_ind += 1

                            # -real(v)'*imag(Q)*imag(v)
                            rate_constraint_qr[q_ind] = real_ind1
                            rate_constraint_qc[q_ind] = imag_ind2
                            rate_constraint_qv[q_ind] = -get_Gi(mm1, mm2)
                            q_ind += 1

                            # imag(v)'*imag(Q)*real(v)
                            rate_constraint_qr[q_ind] = imag_ind1
                            rate_constraint_qc[q_ind] = real_ind2
                            rate_constraint_qv[q_ind] = get_Gi(mm1, mm2)
                            q_ind += 1
                        end
                    end
                end
            end
        end

        # Constant part
        rhs = logabsdet(state.W[k]) - abs(trace(state.W[k]*(UniformScaling(1.) + sigma2s[k]*state.U[k]'*state.U[k]))) + d

        Gurobi.add_qconstr!(model,
            rate_constraint_lind,
            rate_constraint_lval,
            rate_constraint_qr,
            rate_constraint_qc,
            rate_constraint_qv,
            '<',
            rhs)
    end; end

    # Power constraints
    power_constraint_qr = Array(Int, num_vars_per_power_constraint)
    power_constraint_qc = Array(Int, num_vars_per_power_constraint)
    power_constraint_qv = Array(Float64, num_vars_per_power_constraint)
    for i = 1:I
        q_ind = 1
        for k = served_MS_ids(i, assignment)
            for m = 1:M; for n = 1:d
                real_ind = v_ind(k, m, n, false)
                power_constraint_qr[q_ind] = real_ind
                power_constraint_qc[q_ind] = real_ind
                power_constraint_qv[q_ind] = 1.
                q_ind += 1

                imag_ind = v_ind(k, m, n, true)
                power_constraint_qr[q_ind] = imag_ind
                power_constraint_qc[q_ind] = imag_ind
                power_constraint_qv[q_ind] = 1.
                q_ind += 1
            end; end
        end

        Gurobi.add_qconstr!(model,
            [], # lind
            [], # lval
            power_constraint_qr,
            power_constraint_qc,
            power_constraint_qv,
            '<',
            Ps[i])
    end

    # Solve it!
    Gurobi.update_model!(model)
    Gurobi.optimize(model)

    if Gurobi.get_status(model) == :optimal || Gurobi.get_status(model) == :suboptimal
        if Gurobi.get_status(model) == :suboptimal
            Lumberjack.warn("Suboptimal precoder solution found in Razaviyayn2013_MinMaxWMMSE")
        end

        sol = Gurobi.get_solution(model)

        # Put into state variable
        for i = 1:channel.I
            for k in served_MS_ids(i, assignment)
                for m = 1:M
                    for n = 1:d
                        real_ind = v_ind(k, m, n, false)
                        imag_ind = v_ind(k, m, n, true)
                        state.V[k][m,n] = sol[real_ind] + im*sol[imag_ind]
                    end
                end
            end
        end
    else
        Lumberjack.warn("Numerical problems with Gurobi in Razaviyayn2013_MinMaxWMMSE", @compat Dict(:gurobi_output => string(Gurobi.get_status(model))))
    end
end
