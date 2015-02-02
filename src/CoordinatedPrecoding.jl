#############################################################################
# CoordinatedPrecoding
# Description here
# See http://github.com/rasmusbrandt/CoordinatedPrecoding.jl
#############################################################################

module CoordinatedPrecoding

import Gurobi, PyPlot, Lumberjack

export

# Channels
    Channel, SinglecarrierChannel, MulticarrierChannel,
    get_channel_gains_dB, get_average_channel_gains_dB,

# Networks
    Network,
    get_aux_precoding_params,
    get_no_MSs, get_no_BSs, get_no_MSs_per_cell,
    get_no_antennas, get_no_MS_antennas, get_no_BS_antennas,
    get_transmit_powers, get_receiver_noise_powers, get_no_streams,
    set_transmit_powers!, set_receiver_noise_powers!, set_no_streams!,
    get_user_priority, get_user_priorities,
    set_user_priority, set_user_priorities,
    require_equal_no_MS_antennas, require_equal_no_BS_antennas,
    require_single_antenna_MSs, require_single_antenna_BSs,
    require_equal_no_streams, require_single_stream,
    InterferingBroadcastChannel, setup_interfering_broadcast_channel,
    Triangular3SiteNetwork, setup_triangular3site_network,
    ITU_R_InH_Network, setup_itu_r_inh_network,
    draw_user_drop!, draw_channel,

# Cell assignment
    CellAssignment,
    serving_BS_id, served_MS_ids, served_MS_ids_except_me,
    assign_cells_by_id, assign_cells_by_pathloss,
    assign_cells_by_instantaneous_channels,
    require_equal_no_MSs_per_cell,

# Precoding
    Eigenprecoding,
    Gomadam2008_MaxSINR,
    Komulainen2013_WMMSE,
    Razaviyayn2013_MinMaxWMMSE,
    Shi2011_WMMSE,
    SumMSEMinimization,
    zero_receivers, initial_precoders, unity_MSE_weights,

# Simulation
    simulate_convergence, process_convergence, plot_convergence,
    simulate_SNR, process_SNR, plot_SNR,
    simulate_performance,

# Utilities
    clean_simulation_params_for_jld

##########################################################################
# Channels
include("channel.jl")
# Networks
include("network/network.jl")
# Cell assignment
include("cell_assignment/cell_assignment.jl")
# Precoding algorithms
include("precoding/precoding.jl")
# Simulation
include("simulation/simulation.jl")
# Other useful stuff
include("utils.jl")
##########################################################################

end
