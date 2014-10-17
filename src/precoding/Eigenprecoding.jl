immutable EigenprecodingState
    V::Array{Matrix{Complex128},1}
end

function Eigenprecoding(channel::SinglecarrierChannel, network::Network,
    cell_assignment::CellAssignment, settings=Dict())

    Ps = get_transmit_powers(network)
    sigma2s = get_receiver_noise_powers(network)

    settings = check_and_defaultize_settings(Gomadam2008_MaxSINRState, settings)

    state = EigenprecodingState(Array(Matrix{Complex128}, channel.K))

    for i = 1:channel.I
        served = served_MS_ids(i, cell_assignment)
        Kc = length(served)

        for k in served
            _, sing_vals, r_sing_vecs = svd(channel.H[k,i], thin=true)
            sort!(sing_vals, rev=true) # just in case
            max_d = length(sing_vals)

            # Due to sorting, the strongest channel is first
            noise_over_channel_power = sigma2s[k]./sing_vals

            # Activate all channels
            active_channels = [1:max_d]

            Palloc = zeros(active_channels)
            while true
                # Current water level and power allocation
                mu = (1/length(active_channels))*(Ps[i]/Kc + sum(noise_over_channel_power[active_channels]))
                Palloc = mu - noise_over_channel_power[active_channels]

                # If water level high enough, done. Otherwise, turn off weakest channel
                Palloc[end] > 0 ? break : pop!(active_channels)
            end

            # Final precoder for this user
            state.V[k] = r_sing_vecs*diagm([sqrt(Palloc), zeros(Float64, max_d - length(Palloc))])
        end
    end

    return calculate_rates(state, channel, sigma2s, cell_assignment, settings)
end

function check_and_defaultize_settings(::Type{EigenprecodingState},
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

function calculate_rates(state::EigenprecodingState,
    channel::SinglecarrierChannel, sigma2s::Vector{Float64},
    cell_assignment::CellAssignment, settings)

    max_d = min(maximum(channel.Ns), maximum(channel.Ms)) # might not be tight..

    if settings["output_protocol"] == 1
        intercell_tdma_rates = Array(Float64, channel.K, max_d, settings["stop_crit"])
        intracell_tdma_rates = Array(Float64, channel.K, max_d, settings["stop_crit"])
        uncoord_rates = Array(Float64, channel.K, max_d, settings["stop_crit"])
    elseif settings["output_protocol"] == 2
        intercell_tdma_rates = Array(Float64, channel.K, max_d)
        intracell_tdma_rates = Array(Float64, channel.K, max_d)
        uncoord_rates = Array(Float64, channel.K, max_d)
    end

    for i = 1:channel.I
        served = served_MS_ids(i, cell_assignment)
        Kc = length(served)

        for k in served
            Phi_intracell = Hermitian(complex(sigma2s[k]*eye(channel.Ns[k])))
            for j in delete!(IntSet(1:channel.I), i)
                for l in served_MS_ids(j, cell_assignment)
                    #Phi_intracell += Hermitian(channel.H[k,j]*(state.V[l]*state.V[l]')*channel.H[k,j]')
                    herk!(Phi_intracell.uplo, 'N', complex(1.), channel.H[k,j]*state.V[l], complex(1.), Phi_intracell.S)
                end
            end
            Phi_uncoord = copy(Phi_intracell)
            for l in served_MS_ids_except_me(k, i, cell_assignment)
                #Phi_uncoord += Hermitian(channel.H[k,j]*(state.V[l]*state.V[l]')*channel.H[k,j]')
                herk!(Phi_uncoord.uplo, 'N', complex(1.), channel.H[k,i]*state.V[l], complex(1.), Phi_uncoord.S)
            end

            d = size(state.V[k], 2)
            r_intercell = (1/channel.K)*log2(abs(eigvals(eye(d) + (channel.K/Kc)*state.V[k]'*channel.H[k,i]'*(1/sigma2s[k])*channel.H[k,i]*state.V[k])))
            r_intracell = (1/Kc)*log2(abs(eigvals(eye(d) + state.V[k]'*channel.H[k,i]'*(Phi_intracell\channel.H[k,i])*state.V[k])))
            r_uncoord = log2(abs(eigvals(eye(d) + state.V[k]'*channel.H[k,i]'*(Phi_uncoord\channel.H[k,i])*state.V[k])))

            if settings["output_protocol"] == 1
                for iter = 1:settings["stop_crit"]
                    intercell_tdma_rates[k,:,iter] = cat(1, r_intercell, zeros(Float64, max_d - d))
                    intracell_tdma_rates[k,:,iter] = cat(1, r_intracell, zeros(Float64, max_d - d))
                    uncoord_rates[k,:,iter] = cat(1, r_uncoord, zeros(Float64, max_d - d))
                end
            elseif settings["output_protocol"] == 2
                intercell_tdma_rates[k,:] = cat(1, r_intercell, zeros(Float64, max_d - d))
                intracell_tdma_rates[k,:] = cat(1, r_intracell, zeros(Float64, max_d - d))
                uncoord_rates[k,:] = cat(1, r_uncoord, zeros(Float64, max_d - d))
            end
        end
    end

    return [ "intercell_tdma_rates" => intercell_tdma_rates,
             "intracell_tdma_rates" => intracell_tdma_rates,
             "uncoord_rates" => uncoord_rates ]
end
