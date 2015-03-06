##########################################################################
# The RandomLargeScaleNetwork is a network whose geography changes in
# each user drop. The BSs are uniformly distributed over a rectangle,
# and the MSs are distributed on a circle with a certain distance from
# the serving BS.

##########################################################################
# Network definition
type RandomLargeScaleNetwork{MS_t <: PhysicalMS, BS_t <: PhysicalBS, System_t <: System, PropagationEnvironment_t <: PropagationEnvironment} <: PhysicalNetwork
    MSs::Vector{MS_t}
    BSs::Vector{BS_t}

    system::System_t
    no_MSs_per_cell::Int
    propagation_environment::PropagationEnvironment_t
    geography_width::Float64
    geography_height::Float64
    MS_serving_BS_distance::Float64

    assignment::Assignment
end

# Convenience constructor without assignments
RandomLargeScaleNetwork(MSs, BSs, system, no_MSs_per_cell, propagation_environment, geography_width, geography_height, MS_serving_BS_distance) =
    RandomLargeScaleNetwork(MSs, BSs, system, no_MSs_per_cell, propagation_environment, geography_width, geography_height, MS_serving_BS_distance, Assignment())

Base.show(io::IO, x::RandomLargeScaleNetwork) =
    print(io, "RandomLargeScaleNetwork(I = $(length(x.BSs)), Kc = $(x.no_MSs_per_cell), w = $(x.geography_width), h = $(x.geography_height), d = $(x.MS_serving_BS_distance))")
Base.showcompact(io::IO, x::RandomLargeScaleNetwork) =
    print(io, "RandomLargeScaleNetwork($(length(x.BSs)), $(x.no_MSs_per_cell), $(x.geography_width), $(x.geography_height), $(x.MS_serving_BS_distance))")

function setup_random_large_scale_network(
    no_BSs, no_MSs_per_cell, no_MS_antennas, no_BS_antennas;
    system = SinglecarrierSystem(),
    propagation_environment = SimpleLargescaleFadingEnvironment(37.6, 15.3, 0, 8),
    geography_width = 1000.,
    geography_height = 1000.,
    MS_serving_BS_distance = 50.,
    transmit_power = 10^(18.2/10), transmit_powers = transmit_power*ones(Float64, no_BSs),
    BS_antenna_gain_params = [ OmnidirectionalAntennaParams(0) for idx = 1:no_BSs ],
    user_priority = 1., user_priorities = user_priority*ones(Float64, no_BSs*no_MSs_per_cell),
    no_streams = 1, no_streamss = no_streams*ones(Int, no_BSs*no_MSs_per_cell),
    MS_antenna_gain_dB = 0., MS_antenna_gains_dB = MS_antenna_gain_dB*ones(Float64, no_BSs*no_MSs_per_cell),
    receiver_noise_figure = 9., receiver_noise_figures = receiver_noise_figure*ones(Float64, no_BSs*no_MSs_per_cell))

    isa(no_MS_antennas, Vector) || (no_MS_antennas = no_MS_antennas*ones(Int, no_BSs*no_MSs_per_cell))
    isa(no_BS_antennas, Vector) || (no_BS_antennas = no_BS_antennas*ones(Int, no_BSs))

    BSs = [ PhysicalBS(no_BS_antennas[i], Position(0, 0), transmit_powers[i], BS_antenna_gain_params[i]) for i = 1:no_BSs ]
    MSs = [ PhysicalMS(no_MS_antennas[k], Position(0, 0), Velocity(0, 0), user_priorities[k], no_streamss[k], MS_antenna_gains_dB[k], receiver_noise_figures[k], SimpleLargescaleFadingEnvironmentState(zeros(Float64, no_BSs), falses(no_BSs))) for k = 1:no_BSs*no_MSs_per_cell ]

    RandomLargeScaleNetwork(MSs, BSs, system, no_MSs_per_cell, propagation_environment, geography_width, geography_height, MS_serving_BS_distance)
end

##########################################################################
# Standard cell assignment functions
function IDCellAssignment!(channel, network::RandomLargeScaleNetwork)
    Kc = network.no_MSs_per_cell; I = get_no_BSs(network)
    cell_assignment = Array(Int, I*Kc)

    for i = 1:I
        cell_assignment[(i-1)*Kc+1:i*Kc] = i
    end

    network.assignment = Assignment(cell_assignment, I)
end

# Due to construction, we should actually do cell assignment only based on ID.
LargeScaleFadingCellAssignment!(channel, network::RandomLargeScaleNetwork) =
    IDCellAssignment!(channel, network)

##########################################################################
# Simulation functions
function draw_user_drop!{MS_t <: PhysicalMS, BS_t <: PhysicalBS, System_t <: System}(network::RandomLargeScaleNetwork{MS_t, BS_t, System_t, SimpleLargescaleFadingEnvironment})
    I = get_no_BSs(network); Kc = network.no_MSs_per_cell

    # Drop BSs uniformly at random. Also drop the correspondingly served MSs.
    for i = 1:I
        BS_x = network.geography_width*rand(); BS_y = network.geography_height*rand()
        network.BSs[i].position = Position(BS_x, BS_y)

        for k_idx = 1:Kc
            k = (i-1)*Kc + k_idx

            # Position
            theta = 2*pi*rand()
            dx = network.MS_serving_BS_distance*cos(theta)
            dy = network.MS_serving_BS_distance*sin(theta)
            network.MSs[k].position = Position(BS_x + dx, BS_y + dy)

            # Shadow fading (uncorrelated between BSs)
            network.MSs[k].propagation_environment_state =
                SimpleLargescaleFadingEnvironmentState(network.propagation_environment.shadow_sigma_dB*randn(I), falses(I)) # LoS not used in this model
        end
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
    Mx = network.geography_width; My = network.geography_height
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