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

    BSs = CanonicalBS[ CanonicalBS(no_BS_antennas, transmit_power) for k = 1:no_links ]
    MSs = CanonicalMS[ CanonicalMS(no_MS_antennas, BSs[k], no_streams, receiver_noise_power) for k = 1:no_links ]    

    InterferenceChannel(MSs, BSs, system, alpha)
end

##########################################################################
# Simulation functions
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
