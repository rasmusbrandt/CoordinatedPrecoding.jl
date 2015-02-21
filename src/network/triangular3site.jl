##########################################################################
# The Triangular3SiteNetwork is a macrocell network with three sites,
# placed on the vertices of an equilateral triangle. The simulation
# parameters are based on 3GPP Case 1 (TR 25.814 and TR 36.814) macrocell
# simulation environment. Distance-dependent pathloss and log-normal
# shadow fading is used. The shadow fading is correlated (with correlation
# coefficient 0.5) from each MS to the three BSs, but independent
# between MSs.
#
# The typical transmit power is 46 dBm over 10 MHz (600 subcarriers).
# This leads to 10*log10(10^(46/10)/600) = 18.2 dBm for one 15 kHz subcarrier.

##########################################################################
# Six sector antenna params from 3gpp TR 25.996
immutable SixSector3gppAntennaParams <: AntennaParams
    bore_sight_angle::Float64
    angle_with_3dB_gain::Float64
    min_antenna_gain_dB::Float64
end

get_antenna_gain(antenna_params::SixSector3gppAntennaParams, angle) =
    10^(-min(12*(angle/antenna_params.angle_with_3dB_gain)^2, antenna_params.min_antenna_gain_dB)/10)

get_angle{AntennaParams_t <: SixSector3gppAntennaParams}(MS, BS::PhysicalBS{AntennaParams_t}) =
    get_angle(MS.position, BS.position) - BS.antenna_params.bore_sight_angle

##########################################################################
# Network definition
type Triangular3SiteNetwork{MS_t <: PhysicalMS, BS_t <: PhysicalBS, System_t <: System, PropagationEnvironment_t <: PropagationEnvironment} <: PhysicalNetwork
    MSs::Vector{MS_t}
    BSs::Vector{BS_t}

    system::System_t
    no_MSs_per_cell::Int
    propagation_environment::PropagationEnvironment_t
    inter_site_distance::Float64
    guard_distance::Float64

    assignment::Assignment
end

# Convenience constructor without assignments
Triangular3SiteNetwork(MSs, BSs, system, no_MSs_per_cell, propagation_environment, inter_site_distance, guard_distance) =
    Triangular3SiteNetwork(MSs, BSs, system, no_MSs_per_cell, propagation_environment, inter_site_distance, guard_distance, Assignment())

Base.show(io::IO, x::Triangular3SiteNetwork) =
    print(io, "Triangular3Site(I = $(length(x.BSs)), Kc = $(x.no_MSs_per_cell), ISD = $(x.inter_site_distance), GD = $(x.guard_distance))")
Base.showcompact(io::IO, x::Triangular3SiteNetwork) =
    print(io, "Triangular3Site($(length(x.BSs)), $(x.no_MSs_per_cell), $(x.inter_site_distance), $(x.guard_distance))")

# The default parameter values are taken from 3GPP Case 1
# (TR 25.814 and TR 36.814).
function setup_triangular3site_network(
    no_BSs, no_MSs_per_cell, no_MS_antennas, no_BS_antennas;
    system = SinglecarrierSystem(2e9, 15e3),
    propagation_environment = SimpleLargescaleFadingEnvironment(37.6, 15.3, 20, 8),
    inter_site_distance = 500.,
    guard_distance = 35.,
    transmit_power = 10^(18.2/10), transmit_powers = transmit_power*ones(Float64, no_BSs),
    BS_antenna_gain_params =
        [SixSector3gppAntennaParams(-90/180*pi, 35/180*pi, 23),
         SixSector3gppAntennaParams( 30/180*pi, 35/180*pi, 23),
         SixSector3gppAntennaParams(150/180*pi, 35/180*pi, 23)],
    user_priority = 1., user_priorities = user_priority*ones(Float64, no_BSs*no_MSs_per_cell),
    no_streams = 1, no_streamss = no_streams*ones(Int, no_BSs*no_MSs_per_cell),
    MS_antenna_gain_dB = 0., MS_antenna_gains_dB = MS_antenna_gain_dB*ones(Float64, no_BSs*no_MSs_per_cell),
    receiver_noise_figure = 9., receiver_noise_figures = receiver_noise_figure*ones(Float64, no_BSs*no_MSs_per_cell))

    # Consistency check
    if no_BSs != 3
        Lumberjack.error("Triangular3SiteNetwork only allows for I = 3.")
    end

    if !isa(no_MS_antennas, Vector)
        no_MS_antennas = no_MS_antennas*ones(Int, no_BSs*no_MSs_per_cell)
    end
    if !isa(no_BS_antennas, Vector)
        no_BS_antennas = no_BS_antennas*ones(Int, no_BSs)
    end

    BSs = [
        PhysicalBS(no_BS_antennas[1],
            Position(0, sqrt((inter_site_distance/2)^2 + (inter_site_distance/(2*sqrt(3)))^2)),
            transmit_powers[1], BS_antenna_gain_params[1]),
        PhysicalBS(no_BS_antennas[2],
            Position(-inter_site_distance/2, -inter_site_distance/(2*sqrt(3))),
            transmit_powers[2], BS_antenna_gain_params[2]),
        PhysicalBS(no_BS_antennas[3],
            Position(+inter_site_distance/2, -inter_site_distance/(2*sqrt(3))),
            transmit_powers[3], BS_antenna_gain_params[3])
    ]
    MSs = [ PhysicalMS(no_MS_antennas[k], Position(0, 0), Velocity(0, 0), user_priorities[k], no_streamss[k], MS_antenna_gains_dB[k], receiver_noise_figures[k], SimpleLargescaleFadingEnvironmentState(zeros(Float64, no_BSs), falses(no_BSs))) for k = 1:3*no_MSs_per_cell ]

    Triangular3SiteNetwork(MSs, BSs, system, no_MSs_per_cell, 
        propagation_environment, inter_site_distance, guard_distance)
end

##########################################################################
# Standard cell assignment functions
function assign_cells_by_id!{MS_t <: PhysicalMS, BS_t <: PhysicalBS, System_t <: System, PropagationEnvironment_t <: PropagationEnvironment}(network::Triangular3SiteNetwork{MS_t,BS_t,System_t,PropagationEnvironment_t})
    Kc = network.no_MSs_per_cell; I = get_no_BSs(network)
    cell_assignment = Array(Int, I*Kc)

    for i = 1:I
        cell_assignment[(i-1)*Kc+1:i*Kc] = i
    end

    network.assignment = Assignment(cell_assignment, I)
end

##########################################################################
# Simulation functions
function draw_user_drop!{MS_t <: PhysicalMS, BS_t <: PhysicalBS, System_t <: System}(network::Triangular3SiteNetwork{MS_t, BS_t, System_t, SimpleLargescaleFadingEnvironment})
    I = get_no_BSs(network); K = get_no_MSs(network)

    # Shadow fading covariance (correlation coefficient 0.5 between cells)
    Cov_shadow_sqrtm = chol(network.propagation_environment.shadow_sigma_dB^2*(ones(3,3) + eye(3))/2)

    for k = 1:K
        # Generate user position within standard triangle [0,30] degrees, to
        # the right of the base. Apply guard distance to x coordinate.
        xtri = sqrt(((network.inter_site_distance/2)^2 - network.guard_distance^2)*rand() + network.guard_distance^2)
        ytri = (xtri/sqrt(3))*rand()

        # Flip it over with probability 0.5
        if rand() < 0.5
            pos = [xtri;ytri]
        else
            pos = [xtri;-ytri]

            theta = 60*pi/180
            pos = [cos(theta) -sin(theta);sin(theta) cos(theta)]*pos
        end

        # Generate rotation angle
        i = div(k - 1, network.no_MSs_per_cell) + 1 # serving BS id
        if i == 1
            theta = 240*pi/180
        elseif i == 2
            theta = 0.
        elseif i == 3
            theta = 120*pi/180
        end

        # Get final user position
        position = [cos(theta) -sin(theta);sin(theta) cos(theta)]*pos
        network.MSs[k].position = (Position(position[1], position[2]) 
                                   + network.BSs[i].position)

        # Shadow fading
        network.MSs[k].propagation_environment_state =
            SimpleLargescaleFadingEnvironmentState(Cov_shadow_sqrtm*randn(I), falses(I)) # LoS not used in this model
    end
end

function draw_channel{MS_t <: PhysicalMS, BS_t <: PhysicalBS}(network::Triangular3SiteNetwork{MS_t, BS_t, SinglecarrierSystem, SimpleLargescaleFadingEnvironment})
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
function plot_network_layout(network::Triangular3SiteNetwork)
    I = get_no_BSs(network); K = get_no_MSs(network)

    fig = PyPlot.figure()
    ax = fig[:add_subplot](1, 1, 1)

    # Equilateral triangle, with cell edges
    cm = mean([ [network.BSs[i].position.x, network.BSs[i].position.y] for i = 1:I ])
    for i = 1:I
        pos1 = network.BSs[i].position
        pos2 = network.BSs[mod(i, I) + 1].position
        ax[:plot]([pos1.x, pos2.x], [pos1.y, pos2.y]; color="b", linestyle="-")
        ax[:plot]([mean([pos1.x, pos2.x]), cm[1]], [mean([pos1.y, pos2.y]), cm[2]]; color="b", linestyle="--")
    end

    # BSs
    for i = 1:I
        pos = network.BSs[i].position
        ax[:plot](pos.x, pos.y; marker="x", color="b", markersize=10)
    end

    # MSs
    for k = 1:K
        pos = network.MSs[k].position
        ax[:plot](pos.x, pos.y; marker="o", color="r", markersize=10)
    end

    display(fig)
end
