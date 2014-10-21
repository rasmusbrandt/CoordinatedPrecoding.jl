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
                _, _, Vtmp = svd(channel.H[k,i])
                V[k] = sqrt(Ps[i]/Kc)*Vtmp[:, 1:ds[k]]/vecnorm(Vtmp[:, 1:ds[k]])
            end
        end
    end

    return V
end

# Standard algorithms
include("Eigenprecoding.jl")
include("Gomadam2008_MaxSINR.jl")
include("Komulainen2013_WMMSE.jl")
include("Razaviyayn2013_MaxMinWMMSE.jl")
include("Shi2011_WMMSE.jl")
