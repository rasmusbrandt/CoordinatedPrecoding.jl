##########################################################################
# The RandomLargeScaleNetwork is a network whose geography changes in
# each user drop. Both BSs and MSs are uniformly distributed over
# a rectangle.

##########################################################################
# Network definition
type RandomLargeScaleNetwork{MS_t <: PhysicalMS, BS_t <: PhysicalBS, System_t <: System, PropagationEnvironment_t <: PropagationEnvironment} <: PhysicalNetwork
    MSs::Vector{MS_t}
    BSs::Vector{BS_t}

    system::System_t
    propagation_environment::PropagationEnvironment_t
    geography_size::(Float64, Float64)
    aux_network_params::AuxNetworkParams

    assignment::Assignment
end

# Convenience constructor without network params and assignments
RandomLargeScaleNetwork(MSs, BSs, system, propagation_environment, geography_size) =
    RandomLargeScaleNetwork(MSs, BSs, system, propagation_environment, geography_size, AuxNetworkParams(), Assignment())

Base.show(io::IO, x::RandomLargeScaleNetwork) =
    print(io, "RandomLargeScaleNetwork(I = $(length(x.BSs)), I = $(length(x.MSs)), geography = $(x.geography_size))")
Base.showcompact(io::IO, x::RandomLargeScaleNetwork) =
    print(io, "RandomLargeScaleNetwork($(length(x.BSs)), $(length(x.MSs)), $(x.geography_size))")

function setup_random_large_scale_network(
    no_BSs, no_MSs, no_MS_antennas, no_BS_antennas;
    system = SinglecarrierSystem(2e9, 15e3),
    propagation_environment = SimpleLargescaleFadingEnvironment(37.6, 15.3, 0, 8),
    geography_size = (500., 500.),
    transmit_power = 10^(18.2/10), transmit_powers = transmit_power*ones(Float64, no_BSs),
    BS_antenna_gain_params = [ OmnidirectionalAntennaParams(0) for idx = 1:no_BSs ],
    user_priority = 1., user_priorities = user_priority*ones(Float64, no_MSs),
    no_streams = 1, no_streamss = no_streams*ones(Int, no_MSs),
    MS_antenna_gain_dB = 0., MS_antenna_gains_dB = MS_antenna_gain_dB*ones(Float64, no_MSs),
    receiver_noise_figure = 9., receiver_noise_figures = receiver_noise_figure*ones(Float64, no_MSs))

    isa(no_MS_antennas, Vector) || (no_MS_antennas = no_MS_antennas*ones(Int, no_MSs))
    isa(no_BS_antennas, Vector) || (no_BS_antennas = no_BS_antennas*ones(Int, no_BSs))

    BSs = [ PhysicalBS(no_BS_antennas[i], Position(0, 0), transmit_powers[i], BS_antenna_gain_params[i]) for i = 1:no_BSs ]
    MSs = [ PhysicalMS(no_MS_antennas[k], Position(0, 0), Velocity(0, 0), user_priorities[k], no_streamss[k], MS_antenna_gains_dB[k], receiver_noise_figures[k], SimpleLargescaleFadingEnvironmentState(zeros(Float64, no_BSs), falses(no_BSs))) for k = 1:no_MSs ]

    RandomLargeScaleNetwork(MSs, BSs, system, propagation_environment, geography_size)
end

# Greedy scheduler based on the large scale fading realizations
function LargeScaleFadingCellAssignment!(channel, network::RandomLargeScaleNetwork)
    I = get_no_BSs(network); K = get_no_MSs(network)

    aux_params = get_aux_assignment_params(network)
    @defaultize_param! aux_params "max_MSs_per_BS" 1

    # Scheduling matrix
    cell_assignment_matrix = zeros(Int, K, I)

    # User selection metric
    F = (channel.large_scale_fading_factor.^2)*Diagonal(get_transmit_powers(network))
    Fsize = size(F)

    # Do greedy scheduling
    while !all(F .== 0.)
        _, idx = findmax(F)
        k, l = ind2sub(Fsize, idx)

        if sum(cell_assignment_matrix[:,l]) < aux_params["max_MSs_per_BS"]
            cell_assignment_matrix[k,l] = 1
            F[k,:] = 0.
        else
            F[:,l] = 0.
        end
    end

    network.assignment = Assignment(cell_assignment_matrix)

    return AssignmentResults()
end

##########################################################################
# Simulation functions
function draw_user_drop!{MS_t <: PhysicalMS, BS_t <: PhysicalBS, System_t <: System}(network::RandomLargeScaleNetwork{MS_t, BS_t, System_t, SimpleLargescaleFadingEnvironment})
    I = get_no_BSs(network); K = get_no_MSs(network)

    # Drop BSs uniformly at random
    for i = 1:I
        BS_x = network.geography_size[1]*rand(); BS_y = network.geography_size[2]*rand()
        network.BSs[i].position = Position(BS_x, BS_y)
    end

    # Drop MSs uniformly at random
    for k = 1:K
        MS_x = network.geography_size[1]*rand(); MS_y = network.geography_size[2]*rand()
        network.MSs[k].position = Position(MS_x, MS_y)
    end
end

function draw_channel{MS_t <: PhysicalMS, BS_t <: PhysicalBS}(network::RandomLargeScaleNetwork{MS_t, BS_t, SinglecarrierSystem, SimpleLargescaleFadingEnvironment})
    K = get_no_MSs(network); I = get_no_BSs(network)
    Ns = get_no_MS_antennas(network); Ms = get_no_BS_antennas(network)

    coefs = Array(Matrix{Complex128}, K, I)
    large_scale_fading_factor = Array(Float64, K, I)

    distances = get_distances(network)
    angles = get_angles(network)

    for k = 1:K
        for i = 1:I
            # Small scale fading
            coefs[k,i] = (1/sqrt(2))*(randn(Ns[k], Ms[i])
                                 + im*randn(Ns[k], Ms[i]))

            # Pathloss
            pathloss_factor = sqrt(10^(-(network.propagation_environment.pathloss_beta + network.propagation_environment.pathloss_alpha*log10(distances[k,i]))/10))

            # Shadow fading
            shadow_factor = sqrt(10^(network.MSs[k].propagation_environment_state.shadow_realizations_dB[i]/10))

            # Penetration loss
            penetration_loss_factor = sqrt(10^(-network.propagation_environment.penetration_loss_dB/10))

            # BS antenna gain
            bs_antenna_gain = sqrt(get_antenna_gain(network.BSs[i].antenna_params, angles[k,i]))

            # MS antenna gain
            ms_antenna_gain = sqrt(10^(network.MSs[k].antenna_gain_dB/10))

            # Apply scale factors
            large_scale_fading_factor[k,i] = pathloss_factor*shadow_factor*penetration_loss_factor*bs_antenna_gain*ms_antenna_gain
            coefs[k,i] *= large_scale_fading_factor[k,i]
        end
    end

    return SinglecarrierChannel(coefs, Ns, Ms, K, I, large_scale_fading_factor)
end

##########################################################################
# Visualization functions
function plot_network_layout(network::RandomLargeScaleNetwork)
    I = get_no_BSs(network); K = get_no_MSs(network)

    fig = PyPlot.figure()
    ax = fig[:add_subplot](1, 1, 1)

    # Rectangle edges
    Mx = network.geography_size[1]; My = network.geography_size[2]
    ax[:plot]([0, Mx], [0, 0], color="b", linestyle="-")
    ax[:plot]([0, Mx], [My, My], color="b", linestyle="-")
    ax[:plot]([0, 0], [0, My], color="b", linestyle="-")
    ax[:plot]([Mx, Mx], [0, My], color="b", linestyle="-")

    # MSs
    for i = 1:I
        pos = network.BSs[i].position
        ax[:plot](pos.x, pos.y; marker="x", color="b", markersize=6)
        ax[:text](pos.x, pos.y, "BS $i")
    end

    # MSs
    for k = 1:K
        pos = network.MSs[k].position
        ax[:plot](pos.x, pos.y; marker="o", color="r", markersize=6)
        ax[:text](pos.x, pos.y, "MS $k")
    end

    display(fig)
end