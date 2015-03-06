##########################################################################
# Network definition
type InterferingBroadcastChannel{System_t <: System} <: CanonicalNetwork
    MSs::Vector{CanonicalMS}
    BSs::Vector{CanonicalBS}

    system::System_t
    no_MSs_per_cell::Int
    alpha::Float64

    assignment::Assignment
end

# Convenience constructor without assignments
InterferingBroadcastChannel(MSs, BSs, system, no_MSs_per_cell, alpha) =
    InterferingBroadcastChannel(MSs, BSs, system, no_MSs_per_cell, alpha, Assignment())

Base.show(io::IO, x::InterferingBroadcastChannel) =
    print(io, "IBC(I = $(length(x.BSs)), Kc = $(x.no_MSs_per_cell), α = $(x.alpha))")
Base.showcompact(io::IO, x::InterferingBroadcastChannel) =
    print(io, "IBC($(length(x.BSs)), $(x.no_MSs_per_cell), $(x.alpha))")

function setup_interfering_broadcast_channel(
    no_BSs, no_MSs_per_cell, no_MS_antennas, no_BS_antennas;
    system = SinglecarrierSystem(),
    alpha = 1.,
    transmit_power = 1., transmit_powers = transmit_power*ones(Float64, no_BSs),
    user_priority = 1., user_priorities = user_priority*ones(Float64, no_BSs*no_MSs_per_cell),
    no_streams = 1, no_streamss = no_streams*ones(Int, no_BSs*no_MSs_per_cell),
    receiver_noise_power = 1., receiver_noise_powers = receiver_noise_power*ones(Float64, no_BSs*no_MSs_per_cell))

    if !isa(no_MS_antennas, Vector)
        no_MS_antennas = no_MS_antennas*ones(Int, no_BSs*no_MSs_per_cell)
    end
    if !isa(no_BS_antennas, Vector)
        no_BS_antennas = no_BS_antennas*ones(Int, no_BSs)
    end

    BSs = [ CanonicalBS(no_BS_antennas[i], transmit_powers[i]) for i = 1:no_BSs ]
    MSs = [ CanonicalMS(no_MS_antennas[k], user_priorities[k], no_streamss[k], receiver_noise_powers[k]) for k = 1:no_BSs*no_MSs_per_cell ]

    InterferingBroadcastChannel(MSs, BSs, system, no_MSs_per_cell, alpha)
end

##########################################################################
# Standard cell assignment functions
function IDCellAssignment!(channel, network::InterferingBroadcastChannel)
    Kc = network.no_MSs_per_cell; I = get_no_BSs(network)
    cell_assignment = Array(Int, I*Kc)

    for i = 1:I
        cell_assignment[(i-1)*Kc+1:i*Kc] = i
    end

    network.assignment = Assignment(cell_assignment, I)
end

# We don't have large scale fading in this network.
LargeScaleFadingCellAssignment!(channel, network::InterferingBroadcastChannel) =
    IDCellAssignment!(channel, network)

##########################################################################
# Simulation functions
draw_user_drop!(network::InterferingBroadcastChannel) = nothing

function draw_channel{System_t <: SinglecarrierSystem}(network::InterferingBroadcastChannel{System_t})
    I = get_no_BSs(network); Kc = network.no_MSs_per_cell
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
    I = get_no_BSs(network); Kc = network.no_MSs_per_cell; Lc = system.no_subcarriers
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
