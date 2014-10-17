immutable Gomadam2008_MaxSINRState
    U::Array{Matrix{Complex128},1}
    V::Array{Matrix{Complex128},1}

    Phi::Array{Hermitian{Complex128}, 1}
    Gamma::Array{Hermitian{Complex128}, 1}
end

function Gomadam2008_MaxSINR(channel::SinglecarrierChannel, network::Network,
    cell_assignment::CellAssignment, settings=Dict())

    Ps = get_transmit_powers(network)
    sigma2s = get_receiver_noise_powers(network)
    ds = get_no_streams(network)

    settings = check_and_defaultize_settings(Gomadam2008_MaxSINRState, settings)

    state = Gomadam2008_MaxSINRState(
        zero_receivers(channel, ds, cell_assignment),
        initial_precoders(channel, Ps, sigma2s, ds, cell_assignment, settings), 
        Array(Hermitian{Complex128}, channel.K), 
        Array(Hermitian{Complex128}, channel.K))
    rates = Array(Float64, channel.K, maximum(ds), settings["stop_crit"])

    for iter = 1:(settings["stop_crit"]-1)
        update_MSs!(state, channel, sigma2s, ds, cell_assignment)
        rates[:,:,iter] = calculate_rates(channel, state,
                                                  cell_assignment)
        update_BSs!(state, channel, Ps, sigma2s, ds, cell_assignment, settings)
    end
    update_MSs!(state, channel, sigma2s, ds, cell_assignment)
    rates[:,:,end] = calculate_rates(channel, state, cell_assignment)

    if settings["output_protocol"] == 1
        return [ "rates" => rates ]
    elseif settings["output_protocol"] == 2
        return [ "rates" => rates[:,:,end] ]
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

    for k = 1:channel.K
        state.Phi[k] = Hermitian(complex(sigma2s[k]*eye(channel.Ns[k])))
        for j = 1:channel.I
            for l in served_MS_ids(j, cell_assignment)
                #state.Phi[k] += Hermitian(channel.H[k,j]*(state.V[l]*state.V[l]')*channel.H[k,j]')
                herk!(state.Phi[k].uplo, 'N', complex(1.), channel.H[k,j]*state.V[l], complex(1.), state.Phi[k].S)
            end
        end

        # Per-stream receivers
        i = serving_BS_id(k, cell_assignment)
        for d_ind = 1:ds[k]
            Phi_i_plus_n = Hermitian(
                herk!(state.Phi[k].uplo, 'N', complex(-1.), channel.H[k,i]*state.V[k][:,d_ind], complex(1.), copy(state.Phi[k].S)),
                state.Phi[k].uplo)
            ak = Phi_i_plus_n\channel.H[k,i]*state.V[k][:,d_ind]
            state.U[k][:,d_ind] = ak/norm(ak,2)
        end
    end
end

function update_BSs!(state::Gomadam2008_MaxSINRState, channel::SinglecarrierChannel,
    Ps::Vector{Float64}, sigma2s::Vector{Float64}, ds::Vector{Int},
    cell_assignment::CellAssignment, settings)

    for i = 1:channel.I
        state.Gamma[i] = Hermitian(complex(zeros(channel.Ms[i],channel.Ms[i])))
        for k = 1:channel.K
            #state.Gamma[i] += Hermitian(channel.H[k,i]'*(state.U[k]*state.U[k]')*channel.H[k,i])
            herk!(state.Gamma[i].uplo, 'N', complex(1.), channel.H[k,i]'*state.U[k], complex(1.), state.Gamma[i].S)
        end

        # Per-stream precoders
        served = served_MS_ids(i, cell_assignment)
        Nserved = length(served)
        for k in served
            for d_ind = 1:ds[k]
                Gamma_i_plus_n = Hermitian(
                    herk!(state.Gamma[i].uplo, 'N', complex(-1.), channel.H[k,i]'*state.U[k][:,d_ind], complex(1.), copy(state.Gamma[i].S)),
                    state.Gamma[i].uplo)
                vk = (Gamma_i_plus_n + (sigma2s[k]/Ps[i])*eye(channel.Ms[i]))\channel.H[k,i]'*state.U[k][:,d_ind]
                state.V[k][:,d_ind] = sqrt(Ps[i]/(Nserved*ds[k]))*vk/norm(vk,2)
            end
        end
    end
end

function calculate_rates(channel::SinglecarrierChannel,
    state::Gomadam2008_MaxSINRState, cell_assignment::CellAssignment)

    K = length(state.V)
    ds = Int[ size(state.V[k], 2) for k = 1:K ]; max_d = maximum(ds)
    
    rates = Array(Float64, channel.K, max_d)

    for i in 1:channel.I
        for k in served_MS_ids(i, cell_assignment)
            for n = 1:ds[k]
                Fkn = channel.H[k,i]*state.V[k][:,n]
                Phi_i_plus_n = Hermitian(
                    herk!(state.Phi[k].uplo, 'N', complex(-1.), Fkn, complex(1.), copy(state.Phi[k].S)),
                    state.Phi[k].uplo)

                # get rid of imaginary noise
                SINR_term = abs(Fkn'*(Phi_i_plus_n\Fkn))[1]
                rates[k,n] = log2(1 + SINR_term)
            end
        end
    end
    
    return rates
end
