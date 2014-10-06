immutable EigenprecodingState
    V::Array{Matrix{Complex128},1}
end

function Eigenprecoding(channel::SinglecarrierChannel, network::Network,
    cell_assignment::CellAssignment, settings=Dict())

    Ps = get_transmit_powers(network)
    sigma2s = get_receiver_noise_powers(network)

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

    calculate_user_rates(state, channel, sigma2s, cell_assignment)
end

function calculate_user_rates(state::EigenprecodingState,
    channel::SinglecarrierChannel, sigma2s::Vector{Float64},
    cell_assignment::CellAssignment)

    max_d = min(maximum(channel.Ns), maximum(channel.Ms))
    user_rates_intercell = Array(Float64, channel.K, max_d)
    user_rates_intracell = Array(Float64, channel.K, max_d)
    user_rates_uncoord = Array(Float64, channel.K, max_d)

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

            if d < max_d
                user_rates_intercell[k,:] = cat(1, r_intercell, zeros(Float64, max_d))
                user_rates_intracell[k,:] = cat(1, r_intracell, zeros(Float64, max_d))
                user_rates_uncoord[k,:] = cat(1, r_uncoord, zeros(Float64, max_d))
            else
                user_rates_intercell[k,:] = r_intercell
                user_rates_intracell[k,:] = r_intracell
                user_rates_uncoord[k,:] = r_uncoord
            end
        end
    end

    (user_rates_uncoord, user_rates_intercell, user_rates_intracell)
end
