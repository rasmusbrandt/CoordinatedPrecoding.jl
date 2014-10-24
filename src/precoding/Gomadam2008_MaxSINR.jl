immutable Gomadam2008_MaxSINRState
    U::Array{Matrix{Complex128},1}
    W::Array{Hermitian{Complex128},1} # these are only used for rate calculations
    V::Array{Matrix{Complex128},1}
end

function Gomadam2008_MaxSINR(channel::SinglecarrierChannel, network::Network,
    cell_assignment::CellAssignment, settings=Dict())

    settings = check_and_defaultize_settings(Gomadam2008_MaxSINRState, settings)

    Ps = get_transmit_powers(network)
    sigma2s = get_receiver_noise_powers(network)
    ds = get_no_streams(network)

    state = Gomadam2008_MaxSINRState(
        zero_receivers(channel, ds), # needed since I index directly into the matrices right away
        Array(Hermitian{Complex128}, channel.K),
        initial_precoders(channel, Ps, sigma2s, ds, cell_assignment, settings))
    logdet_rates = Array(Float64, channel.K, maximum(ds), settings["stop_crit"])
    MMSE_rates = Array(Float64, channel.K, maximum(ds), settings["stop_crit"])

    for iter = 1:(settings["stop_crit"]-1)
        update_MSs!(state, channel, sigma2s, ds, cell_assignment)
        logdet_rates[:,:,iter] = calculate_logdet_rates(state)
        MMSE_rates[:,:,iter] = calculate_MMSE_rates(state)
        update_BSs!(state, channel, Ps, sigma2s, ds, cell_assignment, settings)
    end
    update_MSs!(state, channel, sigma2s, ds, cell_assignment)
    logdet_rates[:,:,end] = calculate_logdet_rates(state)
    MMSE_rates[:,:,end] = calculate_MMSE_rates(state)

    if settings["output_protocol"] == 1
        return [ "logdet_rates" => logdet_rates, "MMSE_rates" => MMSE_rates ]
    elseif settings["output_protocol"] == 2
        return [ "logdet_rates" => logdet_rates[:,:,end],
                 "MMSE_rates" => MMSE_rates[:,:,end] ]
    end
end

function check_and_defaultize_settings(::Type{Gomadam2008_MaxSINRState},
    settings)

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

    return settings
end

function update_MSs!(state::Gomadam2008_MaxSINRState, channel::SinglecarrierChannel,
    sigma2s::Vector{Float64}, ds::Vector{Int}, cell_assignment::CellAssignment)

    for i = 1:channel.I
        for k in served_MS_ids(i, cell_assignment)
            Phi = Hermitian(complex(sigma2s[k]*eye(channel.Ns[k])))
            for j = 1:channel.I
                for l in served_MS_ids(j, cell_assignment)
                    #Phi += Hermitian(channel.H[k,j]*(state.V[l]*state.V[l]')*channel.H[k,j]')
                    herk!(Phi.uplo, 'N', complex(1.), channel.H[k,j]*state.V[l], complex(1.), Phi.S)
                end
            end

            # Per-stream receivers
            for n = 1:ds[k]
                Phi_i_plus_n = Hermitian(
                    herk!(Phi.uplo, 'N', complex(-1.), channel.H[k,i]*state.V[k][:,n], complex(1.), copy(Phi.S)),
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

function update_BSs!(state::Gomadam2008_MaxSINRState, channel::SinglecarrierChannel,
    Ps::Vector{Float64}, sigma2s::Vector{Float64}, ds::Vector{Int},
    cell_assignment::CellAssignment, settings)

    for i = 1:channel.I
        # Virtual uplink covariance
        Gamma = Hermitian(complex(zeros(channel.Ms[i],channel.Ms[i])))
        for j = 1:channel.I
            for l = served_MS_ids(j, cell_assignment)
                #Gamma += Hermitian(channel.H[k,i]'*(state.U[k]*state.U[k]')*channel.H[k,i])
                herk!(Gamma.uplo, 'N', complex(1.), channel.H[l,i]'*state.U[l], complex(1.), Gamma.S)
            end
        end

        # Per-stream precoders
        served = served_MS_ids(i, cell_assignment)
        Nserved = length(served)
        for k in served
            for n = 1:ds[k]
                Gamma_i_plus_n = Hermitian(
                    herk!(Gamma.uplo, 'N', complex(-1.), channel.H[k,i]'*state.U[k][:,n], complex(1.), copy(Gamma.S)) + (sigma2s[k]/Ps[i])*eye(channel.Ms[i]),
                    Gamma.uplo)
                v = Gamma_i_plus_n\channel.H[k,i]'*state.U[k][:,n]
                state.V[k][:,n] = sqrt(Ps[i]/(Nserved*ds[k]))*v/norm(v,2)
            end
        end
    end
end
