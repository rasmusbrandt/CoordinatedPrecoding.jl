##########################################################################
# The TriangularHetNetNetwork is a cell network with three sites,
# placed on the vertices of an equilateral triangle. The simulation
# parameters are based on 3GPP Case 1 (TR 25.814 and TR 36.814) cell
# simulation environment. Distance-dependent pathloss and log-normal
# shadow fading is used. The shadow fading is correlated (with correlation
# coefficient 0.5) from each MS to the three BSs, but independent
# between MSs.
#
# The typical transmit power is 46 dBm over 10 MHz (600 subcarriers).
# This leads to 10*log10(10^(46/10)/600) = 18.2 dBm for one 15 kHz subcarrier.

##########################################################################
# Network definition
type TriangularHetNetNetwork{MS_t <: PhysicalMS, BS_t <: PhysicalBS, System_t <: System, PropagationEnvironment_t <: PropagationEnvironment} <: PhysicalNetwork
    MSs::Vector{MS_t}
    BSs::Vector{BS_t}

    system::System_t
    no_picos_per_cell::Int
    no_MSs_per_cell::Int
    propagation_environment::PropagationEnvironment_t
    inter_site_distance::Float64
    pico_centre_distance::Float64
    guard_distance::Float64

    assignment::Assignment
end

# Convenience constructor without assignments
TriangularHetNetNetwork(MSs, BSs, system, no_picos_per_cell, no_MSs_per_cell, propagation_environment, inter_site_distance, pico_centre_distance, guard_distance) =
    TriangularHetNetNetwork(MSs, BSs, system, no_picos_per_cell, no_MSs_per_cell, propagation_environment, inter_site_distance, pico_centre_distance, guard_distance, Assignment())

Base.show(io::IO, x::TriangularHetNetNetwork) =
    print(io, "TriangularHetNet(I = $(length(x.BSs)), no_MSs_per_cell = $(x.no_MSs_per_cell), ISD = $(x.inter_site_distance), GD = $(x.guard_distance))")
Base.showcompact(io::IO, x::TriangularHetNetNetwork) =
    print(io, "TriangularHetNet($(length(x.BSs)), $(x.no_MSs_per_cell), $(x.inter_site_distance), $(x.guard_distance))")

# The default parameter values are taken from 3GPP Case 1
# (TR 25.814 and TR 36.814).
function setup_triangularhetnet_network(
    no_picos_per_cell, no_MSs_per_cell, no_MS_antennas, no_BS_antennas;
    system = SinglecarrierSystem(2e9, 15e3),
    propagation_environment = SimpleLargescaleFadingEnvironment(37.6, 15.3, 20, 8),
    inter_site_distance = 500.,
    pico_centre_distance = 100.,
    guard_distance = 35.,
    transmit_power = 10^(18.2/10), transmit_powers = transmit_power*ones(Float64, 3 + 3*no_picos_per_cell),
    BS_antenna_gain_params =
        vcat([SixSector3gppAntennaParams(-deg2rad(90),  deg2rad(35), 23),
              SixSector3gppAntennaParams( deg2rad(30),  deg2rad(35), 23),
              SixSector3gppAntennaParams( deg2rad(150), deg2rad(35), 23)],
             [ OmnidirectionalAntennaParams(0) for idx = 1:3*no_picos_per_cell ]),
    user_priority = 1., user_priorities = user_priority*ones(Float64, 3*no_MSs_per_cell),
    no_streams = 1, no_streamss = no_streams*ones(Int, 3*no_MSs_per_cell),
    MS_antenna_gain_dB = 0., MS_antenna_gains_dB = MS_antenna_gain_dB*ones(Float64, 3*no_MSs_per_cell),
    receiver_noise_figure = 9., receiver_noise_figures = receiver_noise_figure*ones(Float64, 3*no_MSs_per_cell))

    # Consistency check
    if !isa(no_MS_antennas, Vector)
        no_MS_antennas = no_MS_antennas*ones(Int, 3*no_MSs_per_cell)
    end
    if !isa(no_BS_antennas, Vector)
        no_BS_antennas = no_BS_antennas*ones(Int, 3 + 3*no_picos_per_cell)
    end

    BSs = PhysicalBS[]

    # Place macrocells
    macro_centre_distance = (inter_site_distance/2)/cos(deg2rad(30))
    for i = 1:3
        push!(BSs, PhysicalBS(no_BS_antennas[i], rotate(Position(0, macro_centre_distance), deg2rad((i-1)*120)), transmit_powers[i], BS_antenna_gain_params[i]))
    end

    # Place picocells
    no_picos = 3*no_picos_per_cell; pico_rot_angle = deg2rad(360/no_picos)
    pico_base_angle = deg2rad(-(no_picos_per_cell-1)*360/no_picos/2)
    for i = 4:(3 + no_picos)
        push!(BSs, PhysicalBS(no_BS_antennas[i], rotate(Position(0, pico_centre_distance), pico_base_angle + (i-4)*pico_rot_angle), transmit_powers[i], BS_antenna_gain_params[i]))
    end

    MSs = [ PhysicalMS(no_MS_antennas[k], Position(0, 0), Velocity(0, 0), user_priorities[k], no_streamss[k], MS_antenna_gains_dB[k], receiver_noise_figures[k], SimpleLargescaleFadingEnvironmentState(zeros(Float64, 3 + 3*no_picos_per_cell), falses(3 + 3*no_picos_per_cell))) for k = 1:3*no_MSs_per_cell ]

    TriangularHetNetNetwork(MSs, BSs, system, no_picos_per_cell, no_MSs_per_cell, 
        propagation_environment, inter_site_distance, pico_centre_distance, guard_distance)
end

##########################################################################
# The standard cell assignment function only assigns to macrocells.
# Implement some other assignment function to use the picos.
function assign_cells_by_id!{MS_t <: PhysicalMS, BS_t <: PhysicalBS, System_t <: System, PropagationEnvironment_t <: PropagationEnvironment}(channel, network::TriangularHetNetNetwork{MS_t,BS_t,System_t,PropagationEnvironment_t})
    Kc = network.no_MSs_per_cell
    cell_assignment = Array(Int, 3*Kc)

    for i = 1:3
        cell_assignment[(i-1)*Kc+1:i*Kc] = i
    end

    network.assignment = Assignment(cell_assignment, get_no_BSs(network))
end

function assign_cells_by_large_scale_fading!{MS_t <: PhysicalMS, BS_t <: PhysicalBS, System_t <: System, PropagationEnvironment_t <: PropagationEnvironment}(channel, network::TriangularHetNetNetwork{MS_t,BS_t,System_t,PropagationEnvironment_t})
    no_MSs_per_cell = network.no_MSs_per_cell; no_MSs = 3*no_MSs_per_cell
    no_picos_per_cell = network.no_picos_per_cell; no_BSs = 3 + 3*no_picos_per_cell

    aux_params = get_aux_assignment_params(network)
    @defaultize_param! aux_params "TriangularHetNetAssignment:max_MSs_per_BS" 1

    # Scheduling matrix
    cell_assignment_matrix = zeros(Int, no_MSs, no_BSs)

    # User selection metric
    F = (channel.large_scale_fading_factor.^2)*Diagonal(get_transmit_powers(network))
    Fsize = size(F)

    # Do not schedule users in the wrong cell
    for cell = 1:3
        BS_ids = vcat(cell, 3 .+ [(cell-1)*no_picos_per_cell+1:cell*no_picos_per_cell])
        other_cell_MS_ids = setdiff(1:no_MSs, [(cell-1)*no_MSs_per_cell+1:cell*no_MSs_per_cell])
        F[other_cell_MS_ids, BS_ids] = 0.
    end

    # Do greedy scheduling
    while !all(F .== 0.)
        _, idx = findmax(F)
        k, l = ind2sub(Fsize, idx)

        if sum(cell_assignment_matrix[:,l]) < aux_params["TriangularHetNetAssignment:max_MSs_per_BS"]
            cell_assignment_matrix[k,l] = 1
            F[k,:] = 0.
        else
            F[:,l] = 0.
        end
    end

    network.assignment = Assignment(cell_assignment_matrix)
end

##########################################################################
# Simulation functions
function draw_user_drop!{MS_t <: PhysicalMS, BS_t <: PhysicalBS, System_t <: System}(network::TriangularHetNetNetwork{MS_t, BS_t, System_t, SimpleLargescaleFadingEnvironment})
    I = get_no_BSs(network); K = get_no_MSs(network)

    # Shadow fading covariance (correlation coefficient 0.5 between macro cells, no correlation between pico cells)
    correlations = eye(Float64, I, I)
    correlations[1:3, 1:3] = (ones(3,3) + eye(3))/2
    Cov_shadow_sqrtm = network.propagation_environment.shadow_sigma_dB*chol(correlations)

    for k = 1:K
        # Generate user position within standard triangle [0,30] degrees, to
        # the right of the base. Apply guard distance to x coordinate.
        xtri = sqrt(((network.inter_site_distance/2)^2 - network.guard_distance^2)*rand() + network.guard_distance^2)
        ytri = (xtri/sqrt(3))*rand()

        # Flip it over with probability 0.5
        if rand() < 0.5
            pos = Position(xtri, ytri)
        else
            pos = Position(xtri, -ytri)
            pos = rotate(pos, deg2rad(60))
        end

        # Generate rotation angle
        i = div(k - 1, K/3) + 1 # serving BS id
        if i == 1
            θ = deg2rad(240)
        elseif i == 2
            θ = deg2rad(0)
        elseif i == 3
            θ = deg2rad(120)
        end

        # Get final user position
        network.MSs[k].position = network.BSs[i].position + rotate(pos, θ)

        # Shadow fading
        network.MSs[k].propagation_environment_state =
            SimpleLargescaleFadingEnvironmentState(Cov_shadow_sqrtm*randn(I), falses(I)) # LoS not used in this model
    end
end

function draw_channel{MS_t <: PhysicalMS, BS_t <: PhysicalBS}(network::TriangularHetNetNetwork{MS_t, BS_t, SinglecarrierSystem, SimpleLargescaleFadingEnvironment})
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
function plot_network_layout(network::TriangularHetNetNetwork)
    K = get_no_MSs(network)

    fig = PyPlot.figure()
    ax = fig[:add_subplot](1, 1, 1)

    # Equilateral triangle, with cell edges
    cm = mean([ [network.BSs[i].position.x, network.BSs[i].position.y] for i = 1:3 ])
    for i = 1:3
        pos1 = network.BSs[i].position
        pos2 = network.BSs[mod(i, 3) + 1].position
        ax[:plot]([pos1.x, pos2.x], [pos1.y, pos2.y]; color="b", linestyle="-")
        ax[:plot]([mean([pos1.x, pos2.x]), cm[1]], [mean([pos1.y, pos2.y]), cm[2]]; color="b", linestyle="--")
    end

    # Macro BSs
    for i = 1:3
        pos = network.BSs[i].position
        ax[:plot](pos.x, pos.y; marker="x", color="b", markersize=8)
        ax[:text](pos.x, pos.y, "MacroBS $i")

        # Antenna boresight
        spos = Position(network.inter_site_distance/10, 0); epos = rotate(spos, network.BSs[i].antenna_params.bore_sight_angle)
        ax[:plot](pos.x .+ [0, epos.x], pos.y .+ [0, epos.y], color="r", linestyle=":")
    end

    # Pico BSs
    for i = 4:(3 + 3*network.no_picos_per_cell)
        pos = network.BSs[i].position
        ax[:plot](pos.x, pos.y; marker=".", color="b", markersize=6)
        ax[:text](pos.x, pos.y, "PicoBS $i")
    end

    # MSs
    for k = 1:K
        pos = network.MSs[k].position
        ax[:plot](pos.x, pos.y; marker="o", color="r", markersize=6)
        ax[:text](pos.x, pos.y, "MS $k")
    end

    display(fig)
end
