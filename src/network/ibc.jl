##########################################################################
# Network definition
type InterferingBroadcastChannel{System_t <: System} <: CanonicalNetwork
    MSs::Vector{CanonicalMS}
    BSs::Vector{CanonicalBS}

    system::System_t
    no_MSs_per_cell::Int
    alpha::Float64
end
Base.show(io::IO, x::InterferingBroadcastChannel) =
    print(io, "IBC(I = $(length(x.BSs)), Kc = $(x.no_MSs_per_cell), α = $(x.alpha))")
Base.showcompact(io::IO, x::InterferingBroadcastChannel) =
    print(io, "IBC($(length(x.BSs)), $(x.no_MSs_per_cell), $(x.alpha))")

get_no_MSs_per_cell(network::InterferingBroadcastChannel) = network.no_MSs_per_cell

function setup_interfering_broadcast_channel(
    no_BSs::Int, no_MSs_per_cell::Int, no_MS_antennas::Int, no_BS_antennas::Int;
    system = SinglecarrierSystem(),
    alpha::Float64 = 1.,
    transmit_power::Float64 = 1.,
    user_priorities::Vector{Float64} = ones(Float64, no_BSs*no_MSs_per_cell),
    no_streams::Int = 1,
    receiver_noise_power::Float64 = 1.)

    BSs = [ CanonicalBS(no_BS_antennas, transmit_power) for i = 1:no_BSs ]
    MSs = [ CanonicalMS(no_MS_antennas, user_priorities[k], no_streams, receiver_noise_power) for k = 1:no_BSs*no_MSs_per_cell ]

    InterferingBroadcastChannel(MSs, BSs, system, no_MSs_per_cell, alpha)
end

##########################################################################
# Standard cell assignment functions
function assign_cells_by_id{System_t <: System}(network::InterferingBroadcastChannel{System_t})
    Kc = get_no_MSs_per_cell(network); I = get_no_BSs(network)
    assignment = Array(Int, I*Kc)

    for i = 1:I
        assignment[(i-1)*Kc+1:i*Kc] = i
    end

    return CellAssignment(assignment, I)
end

##########################################################################
# Simulation functions
draw_user_drop!(network::InterferingBroadcastChannel) = nothing

function draw_channel{System_t <: SinglecarrierSystem}(network::InterferingBroadcastChannel{System_t})
    I = get_no_BSs(network); Kc = get_no_MSs_per_cell(network)
    Ns = get_no_MS_antennas(network); Ms = get_no_BS_antennas(network)

    coefs = Array(Matrix{Complex128}, I*Kc, I)
    large_scale_fading_factor = Array(Float64, I*Kc, I)

    for k = 1:I*Kc
        serving_id = div(k - 1, Kc) + 1

        for i = 1:I
            if i == serving_id
                large_scale_fading_factor[k,i] = 1
            else
                large_scale_fading_factor[k,i] = sqrt(network.alpha)
            end

            # Apply scale factors
            coefs[k,i] = large_scale_fading_factor[k,i]*(1/sqrt(2))*(randn(Ns[k], Ms[i]) + im*randn(Ns[k], Ms[i]))
        end
    end

    return SinglecarrierChannel(coefs, Ns, Ms, I*Kc, I, large_scale_fading_factor)
end

function draw_channel{System_t <: MulticarrierSystem}(network::InterferingBroadcastChannel{System_t})
    I = get_no_BSs(network); Kc = get_no_MSs_per_cell(network); Lc = system.no_subcarriers
    Ns = get_no_MS_antennas(network); Ms = get_no_BS_antennas(network)

    coefs = Array(Array{Complex128, 3}, I*Kc, I)
    for k = 1:I*Kc
        serving_id = div(k - 1, Kc) + 1

        for i = 1:I
            for l = 1:Lc
                if i == serving_id
                    scale_factor = 1
                else
                    scale_factor = sqrt(network.alpha)
                end

                coefs[k,i,l] = scale_factor*(1/sqrt(2))*(randn(Ns[k], Ms[i]) + im*randn(Ns[k], Ms[i]))
            end
        end
    end

    return MulticarrierChannel(coefs, Ns, Ms, I*Kc, I, Lc)
end
