function zero_receivers(channel::SinglecarrierChannel, ds::Vector{Int},
    cell_assignment::CellAssignment)

    A = Array(Matrix{Complex128}, channel.K)

    for k = 1:channel.K
        A[k] = zeros(Complex128, channel.Ns[k], ds[k])
    end

    return A
end

function initial_precoders(channel::SinglecarrierChannel, Ps::Vector{Float64},
    sigma2s::Vector{Float64}, ds::Vector{Int}, cell_assignment::CellAssignment,
    settings)

    V = Array(Matrix{Complex128}, channel.K)

    if settings["initial_precoders"] == "dft"
        for i = 1:channel.I
            served = served_MS_ids(i, cell_assignment)
            Kc = length(served)

            for k in served
                V[k] = sqrt(Ps[i]/(channel.Ms[i]*ds[k]*Kc))*fft(eye(channel.Ms[i], ds[k]), 1)
            end
        end
    elseif settings["initial_precoders"] == "white"
        for i = 1:channel.I
            served = served_MS_ids(i, cell_assignment)
            Kc = length(served)

            for k in served
                V[k] = sqrt(Ps[i]/(ds[k]*Kc))*eye(channel.Ms[i], ds[k])
            end
        end
    elseif settings["initial_precoders"] == "zeros"
        for i = 1:channel.I
            served = served_MS_ids(i, cell_assignment)
            Kc = length(served)

            for k in served
                V[k] = zeros(channel.Ms[i], ds[k])
            end
        end
    elseif settings["initial_precoders"] == "eigendirection"
        for i = 1:channel.I
            served = served_MS_ids(i, cell_assignment)
            Kc = length(served)

            for k in served
                Vwf = waterfilling(channel.H[k,i], Ps[i]/Kc, sigma2s[k])
                
                # Pad with zeros if necessary
                dalloc = size(Vwf, 2)
                if dalloc < ds[k]
                    V[k] = cat(2, Vwf, zeros(channel.Ms[i], ds[k] - dalloc))
                else
                    V[k] = Vwf
                end
            end
        end
    end

    return V
end

function waterfilling(H::Matrix{Complex128}, P::Float64, sigma2::Float64)
    N, M = size(H)
    _, DD, VV = svd(H, thin=true)

    # The channel strenghts come out sorted, sinc ethe SVD sorts the singular
    # Values. I.e., noise_pow_over_ch_pow starts with the strongest subchannel
    # and ends with the weakest. This is important for the order in which we
    # deactivate subchannels.
    noise_pow_over_ch_pow = sigma2./(DD.^2)

    # Iteratively deactivate subchannels, based on water level
    while true
        # Find waterlevel
        mu = (1/length(noise_pow_over_ch_pow))*(P + sum(noise_pow_over_ch_pow))

        # Get new power allocations
        Psub = mu - noise_pow_over_ch_pow

        # Is the waterlevel high enough?
        if Psub[end] > 0
            break
        else
            # Need to deactive at least one more channel
            pop!(noise_pow_over_ch_pow)
        end
    end

    # Build precoder
    return V = VV[:,1:length(Psub)]*diagm(sqrt(Psub))
end


# Standard algorithms
include("Eigenprecoding.jl")
include("Gomadam2008_MaxSINR.jl")
include("Komulainen2013_WMMSE.jl")
include("Razaviyayn2013_MaxMinWMMSE.jl")
include("Shi2011_WMMSE.jl")
