##########################################################################
# Network definitions
type InterferenceChannel{System_t <: System} <: CanonicalNetwork
    MSs::Vector{CanonicalMS}
    BSs::Vector{CanonicalBS}

    system::System_t
    alpha::Float64
end

type InterferingBroadcastChannel{System_t <: System} <: CanonicalNetwork
    MSs::Vector{CanonicalMS}
    BSs::Vector{CanonicalBS}

    system::System_t
    no_MSs_per_cell::Int
    alpha::Float64
end

get_no_MSs_per_cell(network::InterferingBroadcastChannel) = network.no_MSs_per_cell

function setup_interference_channel(
    no_links::Int, no_MS_antennas::Int, no_BS_antennas::Int;
    system = SinglecarrierSystem(0, 0),
    alpha::Float64 = 1.,
    transmit_power::Float64 = 1.,
    no_streams::Int = 1,
    receiver_noise_power::Float64 = 1.)

    BSs = [ CanonicalBS(no_BS_antennas, transmit_power) for k = 1:no_links ]
    MSs = [ CanonicalMS(no_MS_antennas, no_streams, receiver_noise_power) for k = 1:no_links ]    

    InterferenceChannel(MSs, BSs, system, alpha)
end

function setup_interfering_broadcast_channel(
    no_BSs::Int, no_MSs_per_cell::Int, no_MS_antennas::Int, no_BS_antennas::Int;
    system = SinglecarrierSystem(0, 0),
    alpha::Float64 = 1.,
    transmit_power::Float64 = 1.,
    no_streams::Int = 1,
    receiver_noise_power::Float64 = 1.)

    BSs = [ CanonicalBS(no_BS_antennas, transmit_power) for i = 1:no_BSs ]
    MSs = [ CanonicalMS(no_MS_antennas, no_streams, receiver_noise_power) for k = 1:no_BSs*no_MSs_per_cell ]    

    InterferingBroadcastChannel(MSs, BSs, system, no_MSs_per_cell, alpha)
end

##########################################################################
# Standard cell assignment functions
assign_cells_by_id{System_t <: System}(network::InterferenceChannel{System_t}) =
    CellAssignment([ k for k = 1:get_no_MSs(network) ])

function assign_cells_by_id{System_t <: System}(network::InterferingBroadcastChannel{System_t})
    Kc = get_no_MSs_per_cell(network); I = get_no_BSs(network)
    assignment = Array(Int, I*Kc)

    for i = 1:I
        assignment[(i-1)*Kc+1:i*Kc] = i
    end

    CellAssignment(assignment)
end

##########################################################################
# Simulation functions
draw_user_drop!(network::InterferenceChannel) = nothing
draw_user_drop!(network::InterferingBroadcastChannel) = nothing

function draw_channel{System_t <: SinglecarrierSystem}(network::InterferenceChannel{System_t})

    K = get_no_MSs(network); I = get_no_BSs(network)
    Ns = Int[ network.MSs[k].no_antennas for k = 1:K ]
    Ms = Int[ network.BSs[i].no_antennas for i = 1:I ]

    coefs = Array(Matrix{Complex128}, K, I)
    for k = 1:K
        for i = 1:I
            if k == i
                scale_factor = 1
            else
                scale_factor = sqrt(network.alpha)
            end

            coefs[k,i] = scale_factor*(1/sqrt(2))*(randn(Ns[k], Ms[i]) + im*randn(Ns[k], Ms[i]))
        end
    end

    SinglecarrierChannel(coefs, Ns, Ms, K, I)
end

function draw_channel{System_t <: MulticarrierSystem}(network::InterferenceChannel{System_t})
    
    K = get_no_MSs(network); I = get_no_BSs(network); Lc = system.no_subcarriers
    Ns = Int[ network.MSs[k].no_antennas for k = 1:K ]
    Ms = Int[ network.BSs[i].no_antennas for i = 1:I ]

    coefs = Array(Array{Complex128, 3}, K, I)
    for k = 1:K
        for i = 1:I
            for l = 1:Lc
                if k == i
                    scale_factor = 1
                else
                    scale_factor = sqrt(network.alpha)
                end

                coefs[k,i,l] = scale_factor*(1/sqrt(2))*(randn(Ns[k], Ms[i]) + im*randn(Ns[k], Ms[i]))
            end
        end
    end

    MulticarrierChannel(coefs, Ns, Ms, K, I, Lc)
end

function draw_channel{System_t <: SinglecarrierSystem}(network::InterferingBroadcastChannel{System_t})

    I = get_no_BSs(network); Kc = get_no_MSs_per_cell(network)
    Ns = Int[ network.MSs[k].no_antennas for k = 1:I*Kc ]
    Ms = Int[ network.BSs[i].no_antennas for i = 1:I ]

    coefs = Array(Matrix{Complex128}, I*Kc, I)
    for k = 1:I*Kc
        serving_id = div(k - 1, Kc) + 1

        for i = 1:I
            if i == serving_id
                scale_factor = 1
            else
                scale_factor = sqrt(network.alpha)
            end

            coefs[k,i] = scale_factor*(1/sqrt(2))*(randn(Ns[k], Ms[i]) + im*randn(Ns[k], Ms[i]))
        end
    end

    SinglecarrierChannel(coefs, Ns, Ms, I*Kc, I)
end

function draw_channel{System_t <: MulticarrierSystem}(network::InterferingBroadcastChannel{System_t})
    
    I = get_no_BSs(network); Kc = get_no_MSs_per_cell(network); Lc = system.no_subcarriers
    Ns = Int[ network.MSs[k].no_antennas for k = 1:I*Kc ]
    Ms = Int[ network.BSs[i].no_antennas for i = 1:I ]

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

    MulticarrierChannel(coefs, Ns, Ms, I*Kc, I, Lc)
end
