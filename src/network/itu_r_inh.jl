##########################################################################
# Network definition
type ITU_R_InH_Network{MS_t <: PhysicalMS, BS_t <: PhysicalBS, System_t <: System, PropagationEnvironment_t <: PropagationEnvironment} <: PhysicalNetwork
    MSs::Vector{MS_t}
    BSs::Vector{BS_t}

    system::System_t
    no_MSs_per_cell::Int
    propagation_environments::Vector{PropagationEnvironment_t} # 2 component vector: LoS/NLoS
    inter_site_distance::Float64
    guard_distance::Float64
end

get_no_MSs(network::ITU_R_InH_Network) = 2*network.no_MSs_per_cell
get_no_BSs(network::ITU_R_InH_Network) = 2
get_no_MSs_per_cell(network::ITU_R_InH_Network) = network.no_MSs_per_cell

function setup_itu_r_inh_network{AntennaParams_t <: AntennaParams}(
    no_BSs::Int, no_MSs_per_cell::Int, no_MS_antennas::Int, no_BS_antennas::Int;
    system = SinglecarrierSystem(3.4, 15e3),
    propagation_environments = [SimpleLargescaleFadingEnvironment(16.9, 32.8 + 20*log10(3.4), 0, 3), SimpleLargescaleFadingEnvironment(43.3, 11.5 + 20*log10(3.4), 0, 4)],
    inter_site_distance::Float64 = 60.,
    guard_distance::Float64 = 3.,
    transmit_power::Float64 = 10^(-9.8/10), # 21 dBm over 20 MHz, 1200 used subcarriers
    BS_antenna_gain_params::Vector{AntennaParams_t} = [OmnidirectionalAntennaParams(0), OmnidirectionalAntennaParams(0)],
    no_streams::Int = 1,
    MS_antenna_gain_dB::Float64 = 0.,
    receiver_noise_figure::Float64 = 7.)

    # Consistency check
    if no_BSs != 2
        error("ITU_R_InH_Network only allows for I = 2.")
    end

    BSs = [
        PhysicalBS(no_BS_antennas, Position(0.5*inter_site_distance, 10), transmit_power, BS_antenna_gain_params[1]),
        PhysicalBS(no_BS_antennas, Position(1.5*inter_site_distance, 10), transmit_power, BS_antenna_gain_params[2]) ]
    MSs = [ PhysicalMS(no_MS_antennas, Position(0, 0), Velocity(0, 0), no_streams, MS_antenna_gain_dB, receiver_noise_figure, SimpleLargescaleFadingEnvironmentState(0., false)) for k = 1:2*no_MSs_per_cell ] 

    ITU_R_InH_Network(MSs, BSs, system, no_MSs_per_cell, 
        propagation_environments, inter_site_distance, guard_distance)
end

##########################################################################
# Standard cell assignment functions
function assign_cells_by_id{MS_t <: PhysicalMS, BS_t <: PhysicalBS, System_t <: System, PropagationEnvironment_t <: PropagationEnvironment}(network::ITU_R_InH_Network{MS_t,BS_t,System_t,PropagationEnvironment_t})
    Kc = get_no_MSs_per_cell(network); I = get_no_BSs(network)
    assignment = Array(Int, I*Kc)
inter_site_distance
    for i = 1:I
        assignment[(i-1)*Kc+1:i*Kc] = i
    end

    CellAssignment(assignment)
end

##########################################################################
# Simulation functions
function draw_user_drop!{MS_t <: PhysicalMS, BS_t <: PhysicalBS, System_t <: System}(network::ITU_R_InH_Network{MS_t, BS_t, System_t, SimpleLargescaleFadingEnvironment})
    # ITU-R M.2135-1, p. 33
    P_LoS = d -> begin
        if d <= 18
            return 1
        elseif d < 37
            return exp(-(d - 18)/27)
        else
            return 0.5
        end
    end

    for i = 1:get_no_BSs(network)
        for k = 1:get_no_MSs_per_cell(network)
            # User position?
            while true
                x = network.inter_site_distance*rand(); y = 20*rand()
                i == 2 ? nothing : (x += network.inter_site_distance; y += network.inter_site_distance)
                candidate_position = Position(x, y)

                if get_distance(network.BSs[i].position, candidate_position) > network.guard_distance
                    network.MSs[k].position = candidate_position
                    break
                end
            end

            # LoS/NLoS ?
            if rand() <= P_LoS(get_distance(network.BSs[i].position, network.MSs[k].position))
                network.MSs[k].propagation_environment_state = SimpleLargescaleFadingEnvironmentState(network.propagation_environments[1].shadow_sigma_dB*randn(), true)
            else
                network.MSs[k].propagation_environment_state = SimpleLargescaleFadingEnvironmentState(network.propagation_environments[2].shadow_sigma_dB*randn(), false)
            end
        end
    end
end

function draw_channel{MS_t <: PhysicalMS, BS_t <: PhysicalBS}(network::ITU_R_InH_Network{MS_t, BS_t, SinglecarrierSystem, SimpleLargescaleFadingEnvironment})
    K = get_no_MSs(network); I = get_no_BSs(network)
    Ns = Int[ network.MSs[k].no_antennas for k = 1:K ]
    Ms = Int[ network.BSs[i].no_antennas for i = 1:I ]

    coefs = Array(Matrix{Complex128}, K, 3)

    distances = get_distances(network)

    for k = 1:K
        for i = 1:2
            # Small scale fading
            coefs[k,i] = (1/sqrt(2))*(randn(Ns[k], Ms[i])
                                 + im*randn(Ns[k], Ms[i]))

            # Pathloss
            if network.MSs[k].propagation_environment_state.LoS
                pathloss_factor = sqrt(10^(-(network.propagation_environments[1].pathloss_beta + network.propagation_environments[1].pathloss_alpha*log10(distances[k,i]))/10))
            else
                pathloss_factor = sqrt(10^(-(network.propagation_environments[2].pathloss_beta + network.propagation_environments[2].pathloss_alpha*log10(distances[k,i]))/10))
            end

            # Shadow fading
            shadow_factor = sqrt(10^(network.MSs[k].propagation_environment_state.shadow_realization_dB/10))

            # BS antenna gain
            bs_antenna_gain = sqrt(get_antenna_gain(network.BSs[i].antenna_params))

            # MS antenna gain
            ms_antenna_gain = sqrt(10^(network.MSs[k].antenna_gain_dB/10))

            # Apply scale factors
            coefs[k,i] *= pathloss_factor*shadow_factor*bs_antenna_gain*ms_antenna_gain
        end
    end

    SinglecarrierChannel(coefs, Ns, Ms, K, I)
end

##########################################################################
# Visualization functions
# function plot_network_layout(network::ITU_R_InH_Network)
#     using PyCall; @pyimport matplotlib.pyplot as plt
# end
