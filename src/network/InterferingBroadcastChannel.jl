##########################################################################
# Network definition
type InterferingBroadcastChannel{System_t <: System} <: CanonicalNetwork
    MSs::Vector{CanonicalMS}
    BSs::Vector{CanonicalBS}

    system::System_t
    num_MSs_per_cell::Int
    alpha::Float64
    aux_network_params::AuxNetworkParams

    assignment::Assignment
end

# Convenience constructor without network params and assignments
InterferingBroadcastChannel(MSs, BSs, system, num_MSs_per_cell, alpha) =
    InterferingBroadcastChannel(MSs, BSs, system, num_MSs_per_cell, alpha, AuxNetworkParams(), Assignment())

Base.show(io::IO, x::InterferingBroadcastChannel) =
    print(io, "IBC(I = $(length(x.BSs)), Kc = $(x.num_MSs_per_cell), α = $(x.alpha))")
Base.showcompact(io::IO, x::InterferingBroadcastChannel) =
    print(io, "IBC($(length(x.BSs)), $(x.num_MSs_per_cell), $(x.alpha))")

function setup_interfering_broadcast_channel(
    num_BSs, num_MSs_per_cell, num_MS_antennas, num_BS_antennas;
    system = SinglecarrierSystem(),
    alpha = 1.,
    transmit_power = 1., transmit_powers = transmit_power*ones(Float64, num_BSs),
    user_priority = 1., user_priorities = user_priority*ones(Float64, num_BSs*num_MSs_per_cell),
    num_streams = 1, num_streamss = num_streams*ones(Int, num_BSs*num_MSs_per_cell),
    receiver_noise_power = 1., receiver_noise_powers = receiver_noise_power*ones(Float64, num_BSs*num_MSs_per_cell))

    isa(num_MS_antennas, Vector) || (num_MS_antennas = num_MS_antennas*ones(Int, num_BSs*num_MSs_per_cell))
    isa(num_BS_antennas, Vector) || (num_BS_antennas = num_BS_antennas*ones(Int, num_BSs))

    BSs = [ CanonicalBS(num_BS_antennas[i], transmit_powers[i]) for i = 1:num_BSs ]
    MSs = [ CanonicalMS(num_MS_antennas[k], user_priorities[k], num_streamss[k], receiver_noise_powers[k]) for k = 1:num_BSs*num_MSs_per_cell ]

    InterferingBroadcastChannel(MSs, BSs, system, num_MSs_per_cell, alpha)
end

##########################################################################
# Standard cell assignment functions
function IDCellAssignment!(channel, network::InterferingBroadcastChannel)
    Kc = network.num_MSs_per_cell; I = get_num_BSs(network)
    cell_assignment = Array(Int, I*Kc)

    for i = 1:I
        cell_assignment[(i-1)*Kc+1:i*Kc] = i
    end

    network.assignment = Assignment(cell_assignment, I)

    return AssignmentResults()
end

# We don't have large scale fading in this network.
LargeScaleFadingCellAssignment!(channel, network::InterferingBroadcastChannel) =
    IDCellAssignment!(channel, network)

##########################################################################
# Simulation functions
draw_user_drop!(network::InterferingBroadcastChannel) = nothing

function draw_channel{System_t <: SinglecarrierSystem}(network::InterferingBroadcastChannel{System_t})
    I = get_num_BSs(network); Kc = network.num_MSs_per_cell
    Ns = get_num_MS_antennas(network); Ms = get_num_BS_antennas(network)

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
    I = get_num_BSs(network); Kc = network.num_MSs_per_cell; Lc = system.num_subcarriers
    Ns = get_num_MS_antennas(network); Ms = get_num_BS_antennas(network)

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
