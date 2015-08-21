immutable EigenprecodingState
    V::Array{Matrix{Complex128},1}
end

function Eigenprecoding(channel::SinglecarrierChannel, network)
    assignment = get_assignment(network)

    I = get_num_BSs(network)
    Ps = get_transmit_powers(network)
    sigma2s = get_receiver_noise_powers(network)
    alphas = get_user_priorities(network)
    aux_params = get_aux_precoding_params(network)

    state = EigenprecodingState(Array(Matrix{Complex128}, channel.K))

    for i = 1:I
        served = served_MS_ids(i, assignment)
        Kc = length(served)

        for k in served
            _, sing_vals, r_sing_vecs = svd(channel.H[k,i], thin=true)
            sort!(sing_vals, rev=true) # just in case
            max_d = length(sing_vals)

            # Due to sorting, the strongest channel is first
            noise_over_channel_power = sigma2s[k]./sing_vals

            # Activate all channels
            active_channels = collect(1:max_d)

            Palloc = zeros(Float64, max_d)
            while true
                # Current water level and power allocation
                mu = (1/length(active_channels))*(Ps[i]/Kc + sum(noise_over_channel_power[active_channels]))
                Palloc = mu - noise_over_channel_power[active_channels]

                # If water level high enough, done. Otherwise, turn off weakest channel
                Palloc[end] > 0 ? break : pop!(active_channels)
            end

            # Final precoder for this user
            state.V[k] = r_sing_vecs*diagm(vcat(sqrt(Palloc), zeros(Float64, max_d - length(Palloc))))
        end
    end

    return calculate_logdet_rates(state, channel, sigma2s, alphas, assignment, aux_params)
end

function calculate_logdet_rates(state::EigenprecodingState,
    channel::SinglecarrierChannel, sigma2s, alphas, assignment, aux_params)

    max_d = min(maximum(channel.Ns), maximum(channel.Ms)) # might not be tight..

    if aux_params["output_protocol"] == :all_iterations
        intercell_tdma_logdet_rates = Array(Float64, channel.K, max_d, aux_params["max_iters"])
        intracell_tdma_logdet_rates = Array(Float64, channel.K, max_d, aux_params["max_iters"])
        uncoord_logdet_rates = Array(Float64, channel.K, max_d, aux_params["max_iters"])

        intercell_tdma_MMSE_rates = Array(Float64, channel.K, max_d, aux_params["max_iters"])
        intracell_tdma_MMSE_rates = Array(Float64, channel.K, max_d, aux_params["max_iters"])
        uncoord_MMSE_rates = Array(Float64, channel.K, max_d, aux_params["max_iters"])

        intercell_tdma_weighted_logdet_rates = Array(Float64, channel.K, max_d, aux_params["max_iters"])
        intracell_tdma_weighted_logdet_rates = Array(Float64, channel.K, max_d, aux_params["max_iters"])
        uncoord_weighted_logdet_rates = Array(Float64, channel.K, max_d, aux_params["max_iters"])

        intercell_tdma_weighted_MMSE_rates = Array(Float64, channel.K, max_d, aux_params["max_iters"])
        intracell_tdma_weighted_MMSE_rates = Array(Float64, channel.K, max_d, aux_params["max_iters"])
        uncoord_weighted_MMSE_rates = Array(Float64, channel.K, max_d, aux_params["max_iters"])
    elseif aux_params["output_protocol"] == :final_iteration
        intercell_tdma_logdet_rates = Array(Float64, channel.K, max_d)
        intracell_tdma_logdet_rates = Array(Float64, channel.K, max_d)
        uncoord_logdet_rates = Array(Float64, channel.K, max_d)

        intercell_tdma_MMSE_rates = Array(Float64, channel.K, max_d)
        intracell_tdma_MMSE_rates = Array(Float64, channel.K, max_d)
        uncoord_MMSE_rates = Array(Float64, channel.K, max_d)

        intercell_tdma_weighted_logdet_rates = Array(Float64, channel.K, max_d)
        intracell_tdma_weighted_logdet_rates = Array(Float64, channel.K, max_d)
        uncoord_weighted_logdet_rates = Array(Float64, channel.K, max_d)

        intercell_tdma_weighted_MMSE_rates = Array(Float64, channel.K, max_d)
        intracell_tdma_weighted_MMSE_rates = Array(Float64, channel.K, max_d)
        uncoord_weighted_MMSE_rates = Array(Float64, channel.K, max_d)
    end

    for i = 1:channel.I
        served = served_MS_ids(i, assignment)
        Kc = length(served)

        for k in served
            Phi_intracell = complex(sigma2s[k]*eye(channel.Ns[k]))
            Phi_uncoord = complex(sigma2s[k]*eye(channel.Ns[k]))
            for j in delete!(IntSet(1:channel.I), i)
                for l in served_MS_ids(j, assignment)
                    F = channel.H[k,j]*state.V[l]
                    FFh = F*F'
                    Phi_intracell += FFh
                    Phi_uncoord += FFh
                end
            end
            for l in served_MS_ids_except_me(k, i, assignment)
                F = channel.H[k,i]*state.V[l]
                Phi_uncoord += F*F'
            end

            d = size(state.V[k], 2)
            W_intercell = UniformScaling(1.) + (channel.K/Kc)*state.V[k]'*channel.H[k,i]'*(1/sigma2s[k])*channel.H[k,i]*state.V[k]
            W_intracell = UniformScaling(1.) + state.V[k]'*channel.H[k,i]'*(Phi_intracell\channel.H[k,i])*state.V[k]
            W_uncoord = UniformScaling(1.) + state.V[k]'*channel.H[k,i]'*(Phi_uncoord\channel.H[k,i])*state.V[k]

            r_intercell_logdet = (1/channel.K)*log2(max(1, abs(eigvals(W_intercell))))
            r_intracell_logdet = (1/Kc)*log2(max(1, abs(eigvals(W_intracell))))
            r_uncoord_logdet = log2(max(1, abs(eigvals(W_uncoord))))

            r_intercell_MMSE = (1/channel.K)*log2(max(1, abs(1./diag(inv(W_intercell)))))
            r_intracell_MMSE = (1/Kc)*log2(max(1, abs(1./diag(inv(W_intracell)))))
            r_uncoord_MMSE = log2(max(1, abs(1./diag(inv(W_uncoord)))))

            r_intercell_weighted_logdet = alphas[k]*r_intercell_logdet
            r_intracell_weighted_logdet = alphas[k]*(1/Kc)*r_intracell_logdet
            r_uncoord_weighted_logdet = alphas[k]*r_uncoord_logdet

            r_intercell_weighted_MMSE = alphas[k]*r_intercell_MMSE
            r_intracell_weighted_MMSE = alphas[k]*r_intracell_MMSE
            r_uncoord_weighted_MMSE = alphas[k]*r_uncoord_MMSE

            if aux_params["output_protocol"] == :all_iterations
                for iter = 1:aux_params["max_iters"]
                    intercell_tdma_logdet_rates[k,:,iter] = cat(1, r_intercell_logdet, zeros(Float64, max_d - d))
                    intracell_tdma_logdet_rates[k,:,iter] = cat(1, r_intracell_logdet, zeros(Float64, max_d - d))
                    uncoord_logdet_rates[k,:,iter] = cat(1, r_uncoord_logdet, zeros(Float64, max_d - d))

                    intercell_tdma_MMSE_rates[k,:,iter] = cat(1, r_intercell_MMSE, zeros(Float64, max_d - d))
                    intracell_tdma_MMSE_rates[k,:,iter] = cat(1, r_intracell_MMSE, zeros(Float64, max_d - d))
                    uncoord_MMSE_rates[k,:,iter] = cat(1, r_uncoord_MMSE, zeros(Float64, max_d - d))

                    intercell_tdma_weighted_logdet_rates[k,:,iter] = cat(1, r_intercell_weighted_logdet, zeros(Float64, max_d - d))
                    intracell_tdma_weighted_logdet_rates[k,:,iter] = cat(1, r_intracell_weighted_logdet, zeros(Float64, max_d - d))
                    uncoord_weighted_logdet_rates[k,:,iter] = cat(1, r_uncoord_weighted_logdet, zeros(Float64, max_d - d))

                    intercell_tdma_weighted_MMSE_rates[k,:,iter] = cat(1, r_intercell_weighted_MMSE, zeros(Float64, max_d - d))
                    intracell_tdma_weighted_MMSE_rates[k,:,iter] = cat(1, r_intracell_weighted_MMSE, zeros(Float64, max_d - d))
                    uncoord_weighted_MMSE_rates[k,:,iter] = cat(1, r_uncoord_weighted_MMSE, zeros(Float64, max_d - d))
                end
            elseif aux_params["output_protocol"] == :final_iteration
                intercell_tdma_logdet_rates[k,:] = cat(1, r_intercell_logdet, zeros(Float64, max_d - d))
                intracell_tdma_logdet_rates[k,:] = cat(1, r_intracell_logdet, zeros(Float64, max_d - d))
                uncoord_logdet_rates[k,:] = cat(1, r_uncoord_logdet, zeros(Float64, max_d - d))

                intercell_tdma_MMSE_rates[k,:] = cat(1, r_intercell_MMSE, zeros(Float64, max_d - d))
                intracell_tdma_MMSE_rates[k,:] = cat(1, r_intracell_MMSE, zeros(Float64, max_d - d))
                uncoord_MMSE_rates[k,:] = cat(1, r_uncoord_MMSE, zeros(Float64, max_d - d))

                intercell_tdma_weighted_logdet_rates[k,:] = cat(1, r_intercell_weighted_logdet, zeros(Float64, max_d - d))
                intracell_tdma_weighted_logdet_rates[k,:] = cat(1, r_intracell_weighted_logdet, zeros(Float64, max_d - d))
                uncoord_weighted_logdet_rates[k,:] = cat(1, r_uncoord_weighted_logdet, zeros(Float64, max_d - d))

                intercell_tdma_weighted_MMSE_rates[k,:] = cat(1, r_intercell_weighted_MMSE, zeros(Float64, max_d - d))
                intracell_tdma_weighted_MMSE_rates[k,:] = cat(1, r_intracell_weighted_MMSE, zeros(Float64, max_d - d))
                uncoord_weighted_MMSE_rates[k,:] = cat(1, r_uncoord_weighted_MMSE, zeros(Float64, max_d - d))
            end
        end
    end

    results = PrecodingResults()
    results["intercell_tdma_logdet_rates"] = intercell_tdma_logdet_rates
    results["intracell_tdma_logdet_rates"] = intracell_tdma_logdet_rates
    results["uncoord_logdet_rates"] = uncoord_logdet_rates
    results["intercell_tdma_MMSE_rates"] = intercell_tdma_MMSE_rates
    results["intracell_tdma_MMSE_rates"] = intracell_tdma_MMSE_rates
    results["uncoord_MMSE_rates"] = uncoord_MMSE_rates
    results["intercell_tdma_weighted_logdet_rates"] = intercell_tdma_weighted_logdet_rates
    results["intracell_tdma_weighted_logdet_rates"] = intracell_tdma_weighted_logdet_rates
    results["uncoord_weighted_logdet_rates"] = uncoord_weighted_logdet_rates
    results["intercell_tdma_weighted_MMSE_rates"] = intercell_tdma_weighted_MMSE_rates
    results["intracell_tdma_weighted_MMSE_rates"] = intracell_tdma_weighted_MMSE_rates
    results["uncoord_weighted_MMSE_rates"] = uncoord_weighted_MMSE_rates
    return results
end
