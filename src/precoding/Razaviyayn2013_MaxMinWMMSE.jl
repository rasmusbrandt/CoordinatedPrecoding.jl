immutable Razaviyayn2013_MaxMinWMMSEState
    U::Array{Matrix{Complex128},1}
    W::Array{Hermitian{Complex128},1}
    V::Array{Matrix{Complex128},1}

    Phi::Array{Hermitian{Complex128},1}
    Gamma::Array{Hermitian{Complex128},1}
end

function Razaviyayn2013_MaxMinWMMSE(channel::SinglecarrierChannel,
    network::Network, cell_assignment::CellAssignment, settings=Dict())

    # The implementation is currently limited in the respects below. This is in
    # order to simplify the Gurobi optimization variable indexing. With equal
    # number of antennas and streams, it is very easy to calculate the variable
    # offsets required. (See functions v_ind and mu_ind.)
    require_equal_no_BS_antennas(network)
    require_equal_no_streams(network)
    require_equal_no_MSs_per_cell(cell_assignment)

    Ps = get_transmit_powers(network)
    sigma2s = get_receiver_noise_powers(network)
    ds = get_no_streams(network)

    defaultize_settings!(Razaviyayn2013_MaxMinWMMSEState, settings)

    state = Razaviyayn2013_MaxMinWMMSEState(
        Array(Matrix{Complex128}, channel.K),
        Array(Hermitian{Complex128}, channel.K),
        initial_precoders(channel, Ps, sigma2s, ds, cell_assignment, settings),
        Array(Hermitian{Complex128}, channel.K),
        Array(Hermitian{Complex128}, channel.K))
    user_rates = Array(Float64, channel.K, maximum(ds), settings["stop_crit"])

    for iter = 1:(settings["stop_crit"]-1)
        update_MSs!(state, channel, sigma2s, ds, cell_assignment)
        user_rates[:,:,iter] = calculate_user_rates(state)
        update_BSs!(state, channel, Ps, sigma2s, cell_assignment, settings)
    end
    update_MSs!(state, channel, sigma2s, ds, cell_assignment)
    user_rates[:,:,end] = calculate_user_rates(state)

    return user_rates
end

function defaultize_settings!(::Type{Razaviyayn2013_MaxMinWMMSEState}, settings)
    if !haskey(settings, "stop_crit")
        settings["stop_crit"] = 20
    end
    if !haskey(settings, "initial_precoders")
        settings["initial_precoders"] = "dft"
    end
    if !haskey(settings, "verbose")
        settings["verbose"] = false
    end
end

function update_MSs!(state::Razaviyayn2013_MaxMinWMMSEState,
    channel::SinglecarrierChannel, sigma2s::Vector{Float64}, ds::Vector{Int},
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
        state.U[k] = state.Phi[k]\Fk
        state.W[k] = Hermitian((eye(ds[k]) - state.U[k]'*Fk)\eye(ds[k]))
    end
end

# Some of the vectors describing the power constraints and rate constraints
# could be prebuilt and reused in the iterations. Probably no major performance
# boost would come from that though, based on experience with the
# JointPrecodingMCSSelection implementation.
function update_BSs!(state::Razaviyayn2013_MaxMinWMMSEState,
    channel::SinglecarrierChannel, Ps::Vector{Float64},
    sigma2s::Vector{Float64}, cell_assignment::CellAssignment, settings)

    # The following is OK, since we have checked this at the top of
    # the outer function.
    K = channel.K; I = channel.I; Kc = int(channel.K/channel.I)
    N = channel.Ns[1]; M = channel.Ms[1]; d = size(state.W[1], 2)

    # Bookkeeping for the Gurobi model
    no_t_vars = 1
    no_v_vars = 2*K*M*d
    no_lin_vars_per_rate_constraint = 2*M*d + 1 # including t variable
    no_quad_vars_per_rate_constraint = K*d*(2*M)^2
    no_vars_per_power_constraint = 2*Kc*M*d

    # Index functions for optimization variables
    t_ind = 1
    v_ind(k, m, n, imag::Bool) = (no_t_vars + (k-1)*M*d*2 + (m-1)*d*2 + (n-1)*2 + 1 + int(imag))

    # Gurobi environment
    env = Gurobi.Env()
    if !settings["verbose"]
        Gurobi.setparam!(env, "OutputFlag", 0)
    end
    model = Gurobi.Model(env, "Razaviyayn2013_MaxMinWMMSE", :maximize)

    # Objective
    Gurobi.add_cvar!(model, 1., 0, Inf) # t

    # Other variables not in objective
    Gurobi.add_cvars!(model, zeros(Float64, no_v_vars))

    # Need to update model before adding constraints
    Gurobi.update_model!(model)

    # Rate constraints
    rate_constraint_lind = Array(Int, no_lin_vars_per_rate_constraint)
    rate_constraint_lval = Array(Float64, no_lin_vars_per_rate_constraint)
    rate_constraint_qr = Array(Int, no_quad_vars_per_rate_constraint)
    rate_constraint_qc = Array(Int, no_quad_vars_per_rate_constraint)
    rate_constraint_qv = Array(Float64, no_quad_vars_per_rate_constraint)
    for i = 1:I; for k = served_MS_ids(i, cell_assignment)
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

            for l = served_MS_ids(j, cell_assignment)
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
        rhs = abs(logdet(state.W[k])) - abs(trace(state.W[k]*(eye(d) + sigma2s[k]*state.U[k]'*state.U[k]))) + d

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
    power_constraint_qr = Array(Int, no_vars_per_power_constraint)
    power_constraint_qc = Array(Int, no_vars_per_power_constraint)
    power_constraint_qv = Array(Float64, no_vars_per_power_constraint)
    for i = 1:I
        q_ind = 1
        for k = served_MS_ids(i, cell_assignment)
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

    if Gurobi.get_status(model) == :optimal
        sol = Gurobi.get_solution(model)

        # Put into state variable
        for k = 1:channel.K
            for m = 1:M
                for n = 1:d
                    real_ind = v_ind(k, m, n, false)
                    imag_ind = v_ind(k, m, n, true)
                    state.V[k][m,n] = sol[real_ind] + im*sol[imag_ind]
                end
            end
        end
    else
        warn("Numerical problems with Gurobi")
    end
end

function calculate_user_rates(state::Razaviyayn2013_MaxMinWMMSEState)
    K = length(state.W)
    ds = Int[ size(state.W[k], 1) for k = 1:K ]; max_d = maximum(ds)

    user_rates = Array(Float64, K, max_d)

    for k = 1:K
        # W is p.d., but numerically we might get tiny negative eigenvalues.
        r = log2(abs(eigvals(state.W[k])))

        if ds[k] < max_d
            user_rates[k,:] = cat(1, r, zeros(Float64, max_d - ds[k]))
        else
            user_rates[k,:] = r
        end
    end

    return user_rates
end
