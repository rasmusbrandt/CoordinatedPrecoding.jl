##########################################################################
# The IndoorsNetwork is an indoors microcell network inspired by the ITU-R
# InH indoors network. The channel model is simplified; it is just using
# uncorrelated Rayleigh fading. The large-scale fading parameters should
# be close to the InH model though. The number of BSs has also been
# generalized compared to InH. The LoS is independent between BSs, for
# each MS. The same holds for the shadow fading. This is motivated by
# Sec. 3.3.1 in Deliverable 1.1.2 of IST-WINNER II.
#
# The typical transmit power is 21 dBm over 20 MHz (1200 subcarriers).
# This leads to 10*log10(10^(21/10)/1200) = -9.8 dBm for one 15 kHz subcarrier.

##########################################################################
# Network definition
type IndoorsNetwork{MS_t <: PhysicalMS, BS_t <: PhysicalBS, System_t <: System, PropagationEnvironment_t <: PropagationEnvironment} <: PhysicalNetwork
    MSs::Vector{MS_t}
    BSs::Vector{BS_t}

    system::System_t
    no_MSs_per_cell::Int
    propagation_environments::Dict{Symbol,PropagationEnvironment_t} # two keys only, :LoS and :NLoS
    corridor_length::Float64
    corridor_width::Float64
    guard_distance::Float64

    assignment::Assignment
end

# Convenience constructor without assignments
IndoorsNetwork(MSs, BSs, system, no_MSs_per_cell, propagation_environments, corridor_length, corridor_width, guard_distance) =
    IndoorsNetwork(MSs, BSs, system, no_MSs_per_cell, propagation_environments, corridor_length, corridor_width, guard_distance, Assignment())

Base.show(io::IO, x::IndoorsNetwork) =
    print(io, "Indoors(I = $(length(x.BSs)), Kc = $(x.no_MSs_per_cell), corridor_length = $(x.corridor_length), corridor_width = $(x.corridor_width), GD = $(x.guard_distance))")
Base.showcompact(io::IO, x::IndoorsNetwork) =
    print(io, "Indoors($(length(x.BSs)), $(x.no_MSs_per_cell), $(x.corridor_length), $(x.corridor_width), $(x.guard_distance))")

# The default parameter values are taken from ITU-R M.2135-1.
function setup_indoors_network(
    no_BSs, no_MSs_per_cell, no_MS_antennas, no_BS_antennas;
    system = SinglecarrierSystem(3.4e9, 15e3),
    propagation_environments = [:LoS => SimpleLargescaleFadingEnvironment(16.9, 32.8 + 20*log10(3.4), 0, 3),
                                :NLoS => SimpleLargescaleFadingEnvironment(43.3, 11.5 + 20*log10(3.4), 0, 4)],
    corridor_length = 120.,
    corridor_width = 50.,
    guard_distance = 3.,
    transmit_power = 10^(-9.8/10), transmit_powers = transmit_power*ones(Float64, no_BSs),
    BS_antenna_gain_params = [ OmnidirectionalAntennaParams(0) for i = 1:no_BSs ],
    user_priority = 1., user_priorities = ones(Float64, no_BSs*no_MSs_per_cell),
    no_streams = 1, no_streamss = no_streams*ones(Int, no_BSs*no_MSs_per_cell),
    MS_antenna_gain_dB = 0., MS_antenna_gains_dB = MS_antenna_gain_dB*ones(Float64, no_BSs*no_MSs_per_cell),
    receiver_noise_figure = 7., receiver_noise_figures = receiver_noise_figure*ones(Float64, no_BSs*no_MSs_per_cell))

    if !isa(no_MS_antennas, Vector)
        no_MS_antennas = no_MS_antennas*ones(Int, no_BSs*no_MSs_per_cell)
    end
    if !isa(no_BS_antennas, Vector)
        no_BS_antennas = no_BS_antennas*ones(Int, no_BSs)
    end

    # Place BSs
    y = corridor_width/2
    Δ = corridor_length/no_BSs
    BSs = Array(PhysicalBS, 0)
    for i = 1:no_BSs
        # BSs uniformly placed along corridor
        push!(BSs, PhysicalBS(no_BS_antennas[i], Position((i-1)*Δ + Δ/2, y), transmit_powers[i], BS_antenna_gain_params[i]))
    end

    MSs = [ PhysicalMS(no_MS_antennas[k], Position(0, 0), Velocity(0, 0), user_priorities[k], no_streamss[k], MS_antenna_gains_dB[k], receiver_noise_figures[k], SimpleLargescaleFadingEnvironmentState(zeros(Float64, no_BSs), falses(no_BSs))) for k = 1:no_BSs*no_MSs_per_cell ]

    IndoorsNetwork(MSs, BSs, system, no_MSs_per_cell, 
        propagation_environments, corridor_length, corridor_width, guard_distance)
end

##########################################################################
# Standard cell assignment functions
function IDCellAssignment!{MS_t <: PhysicalMS, BS_t <: PhysicalBS, System_t <: System, PropagationEnvironment_t <: PropagationEnvironment}(channel, network::IndoorsNetwork{MS_t,BS_t,System_t,PropagationEnvironment_t})
    Kc = network.no_MSs_per_cell; I = get_no_BSs(network)
    cell_assignment = Array(Int, I*Kc)

    for i = 1:I
        cell_assignment[(i-1)*Kc+1:i*Kc] = i
    end

    network.assignment = Assignment(cell_assignment, I)
end

##########################################################################
# Simulation functions
function draw_user_drop!{MS_t <: PhysicalMS, BS_t <: PhysicalBS, System_t <: System}(network::IndoorsNetwork{MS_t, BS_t, System_t, SimpleLargescaleFadingEnvironment})
    I = get_no_BSs(network); K = get_no_MSs(network)
    Δ = network.corridor_length/I

    # ITU-R M.2135-1, p. 33
    P_LoS = d -> begin
        if d <= 18
            return 1.
        elseif d < 37
            return exp(-(d - 18)/27)
        else
            return 0.5
        end
    end

    for k = 1:K
        i = div(k - 1, network.no_MSs_per_cell) + 1 # serving BS id

        # MS position
        while true
            x = (i-1)*Δ + Δ*rand()
            y = network.corridor_width*rand()
            candidate_position = Position(x, y)

            if get_distance(network.BSs[i].position, candidate_position) > network.guard_distance
                network.MSs[k].position = candidate_position
                break
            end
        end

        # Shadow fading and LoS/NLoS. Uncorrelated between BS links.
        shadow_realizations_dB = zeros(Float64, I)
        LoSs = BitArray(I)
        for j = 1:I
            if rand() <= P_LoS(get_distance(network.BSs[j].position, network.MSs[k].position))
                LoSs[j] = true
                shadow_realizations_dB[j] = network.propagation_environments[:LoS].shadow_sigma_dB*randn()
            else
                LoSs[j] = false
                shadow_realizations_dB[j] = network.propagation_environments[:NLoS].shadow_sigma_dB*randn()
            end
        end
        network.MSs[k].propagation_environment_state =
            SimpleLargescaleFadingEnvironmentState(shadow_realizations_dB, LoSs)
    end
end

function draw_channel{MS_t <: PhysicalMS, BS_t <: PhysicalBS}(network::IndoorsNetwork{MS_t, BS_t, SinglecarrierSystem, SimpleLargescaleFadingEnvironment})
    K = get_no_MSs(network); I = get_no_BSs(network)
    Ns = get_no_MS_antennas(network); Ms = get_no_BS_antennas(network)

    coefs = Array(Matrix{Complex128}, K, I)
    large_scale_fading_factor = Array(Float64, K, I)

    distances = get_distances(network)

    for k = 1:K
        for i = 1:I
            # Small scale fading
            coefs[k,i] = (1/sqrt(2))*(randn(Ns[k], Ms[i])
                                 + im*randn(Ns[k], Ms[i]))

            # Pathloss
            if network.MSs[k].propagation_environment_state.LoSs[i]
                pathloss_factor =
                    sqrt(10^(-(network.propagation_environments[:LoS].pathloss_beta + network.propagation_environments[:LoS].pathloss_alpha*log10(distances[k,i]))/10))
            else
                pathloss_factor =
                    sqrt(10^(-(network.propagation_environments[:NLoS].pathloss_beta + network.propagation_environments[:NLoS].pathloss_alpha*log10(distances[k,i]))/10))
            end

            # Shadow fading
            shadow_factor = sqrt(10^(network.MSs[k].propagation_environment_state.shadow_realizations_dB[i]/10))

            # BS antenna gain
            bs_antenna_gain = sqrt(get_antenna_gain(network.BSs[i].antenna_params))

            # MS antenna gain
            ms_antenna_gain = sqrt(10^(network.MSs[k].antenna_gain_dB/10))

            # Apply scale factors
            large_scale_fading_factor[k,i] = pathloss_factor*shadow_factor*bs_antenna_gain*ms_antenna_gain
            coefs[k,i] *= large_scale_fading_factor[k,i]
        end
    end

    return SinglecarrierChannel(coefs, Ns, Ms, K, I, large_scale_fading_factor)
end

##########################################################################
# Visualization functions
function plot_network_layout(network::IndoorsNetwork)
    I = get_no_BSs(network); K = get_no_MSs(network)
    Δ = network.corridor_length/I
    δ = network.corridor_width

    fig = PyPlot.figure()
    ax = fig[:add_subplot](1, 1, 1)

    # Corridor outline
    pos1 = network.BSs[1].position; pos2 = network.BSs[end].position
    sw = (pos1.x - Δ/2, pos1.y - δ/2); se = (pos2.x + Δ/2, pos2.y - δ/2)
    nw = (pos1.x - Δ/2, pos1.y + δ/2); ne = (pos2.x + Δ/2, pos2.y + δ/2)
    ax[:plot]([sw[1], se[1]], [sw[2], se[2]]; color="b", linestyle="-", linewidth=4)
    ax[:plot]([se[1], ne[1]], [se[2], ne[2]]; color="b", linestyle="-", linewidth=4)
    ax[:plot]([ne[1], nw[1]], [ne[2], nw[2]]; color="b", linestyle="-", linewidth=4)
    ax[:plot]([nw[1], sw[1]], [nw[2], sw[2]]; color="b", linestyle="-", linewidth=4)

    # BSs
    for i = 1:I
        pos = network.BSs[i].position
        ax[:plot](pos.x, pos.y; marker="x", color="b", markersize=10)

        # Cell outline
        if i < I
            ax[:plot]([pos.x + Δ/2, pos.x + Δ/2], [pos.y - δ/2, pos.y + δ/2]; color="b", linestyle="--", linewidth=2)
        end
    end

    # MSs
    for k = 1:K
        pos = network.MSs[k].position
        ax[:plot](pos.x, pos.y; marker="o", color="r", markersize=10)
    end

    display(fig)
end
