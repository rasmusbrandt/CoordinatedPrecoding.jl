#############################################################################
# CoordinatedPrecoding
# Description here
# See http://github.com/rasmusbrandt/CoordinatedPrecoding.jl
#############################################################################

module CoordinatedPrecoding

# Necessary dependencies
import PyPlot, Lumberjack, ProgressMeter

# Optional dependencies
if Pkg.installed("Gurobi") != nothing
    import Gurobi
end

export

# Channels
    Channel,

    SinglecarrierChannel, MulticarrierChannel,

    get_channel_gains_dB,
    get_average_channel_gains_dB,

# Networks
    Network,

    InterferingBroadcastChannel, setup_interfering_broadcast_channel,
    IndoorsNetwork, setup_indoors_network,
    RandomLargeScaleNetwork, setup_random_large_scale_network,
    TriangularMacroNetwork, setup_triangularmacro_network,
    TriangularHetNetNetwork, setup_triangularhetnet_network,

    draw_user_drop!, draw_channel, plot_network_layout,

    get_no_MSs, get_no_BSs,

    get_assignment, set_assignment!,
    has_aux_network_param,
    get_aux_network_param, set_aux_network_param!,
    get_aux_network_params, set_aux_network_params!,
    has_aux_precoding_param,
    get_aux_precoding_param, set_aux_precoding_param!,
    get_aux_precoding_params, set_aux_precoding_params!,
    has_aux_assignment_param,
    get_aux_assignment_param, set_aux_assignment_param!,
    get_aux_assignment_params, set_aux_assignment_params!,

    get_no_MS_antennas, get_no_BS_antennas,
    get_transmit_powers, set_transmit_powers!,
    get_transmit_powers_dBm, set_transmit_powers_dBm!,
    get_receiver_noise_powers, set_receiver_noise_powers!,
    get_receiver_noise_powers_dBm, set_receiver_noise_powers_dBm!,

    get_no_streams, set_no_streams!,
    get_user_priorities, set_user_priorities!,

    get_average_SNRs, set_average_SNRs!,
    get_average_SNRs_dB, set_average_SNRs_dB!,

    require_equal_no_MS_antennas,
    require_equal_no_BS_antennas,
    require_single_antenna_MSs,
    require_single_antenna_BSs,
    require_equal_no_streams,
    require_single_stream,

# Assignment
    Assignment,

    IDCellAssignment!,
    LargeScaleFadingCellAssignment!,

    AssignmentResults,
    AuxAssignmentParams,

    active_BSs,
    no_served_MSs,
    serving_BS_id,
    served_MS_ids,
    served_MS_ids_except_me,

    coordinated_BS_ids,
    coordinated_MS_ids,

    require_equal_no_MSs_per_cell,

# Precoding
    Eigenprecoding,
    Gomadam2008_MaxSINR,
    Gomadam2008_MinWLI,
    Komulainen2013_WMMSE,
    Razaviyayn2013_MinMaxWMMSE,
    Shi2011_WMMSE,

    PrecodingResults,
    AuxPrecodingParams,

    initial_receivers,
    initial_MSE_weights,
    initial_precoders,

    calculate_logdet_rates,
    calculate_MMSE_rates,
    calculate_weighted_logdet_rates,
    calculate_weighted_MMSE_rates,
    calculate_allocated_power,

# Simulation
    simulate,
    postprocess,
    plot,

    simulate_assignment_convergence,
    postprocess_assignment_convergence,
    plot_assignment_convergence,

    simulate_precoding_convergence,
    postprocess_precoding_convergence,
    plot_precoding_convergence,

    timing,

# Utilities
    @defaultize_param!,
    clean_simulation_params_for_jld

##########################################################################
# General includes

# Networks
include("network/network.jl")

# Channels
include("channel.jl")

# Simulation
include("simulation/simulation.jl")

# Cell/cluster assignment
include("assignment.jl")

# Precoding algorithms
include("precoding/precoding.jl")

# Miscellaneous utilities
include("utils.jl")

##########################################################################
# Specific includes

# Networks
include("network/InterferingBroadcastChannel.jl")
include("network/IndoorsNetwork.jl")
include("network/RandomLargeScaleNetwork.jl")
include("network/TriangularMacroNetwork.jl")
include("network/TriangularHetNetNetwork.jl")

# Precoding methods
include("precoding/Eigenprecoding.jl")
include("precoding/Gomadam2008_MaxSINR.jl")
include("precoding/Gomadam2008_MinWLI.jl")
include("precoding/Komulainen2013_WMMSE.jl")
include("precoding/Razaviyayn2013_MinMaxWMMSE.jl")
include("precoding/Shi2011_WMMSE.jl")

# Simulation utilities
include("simulation/timing.jl")
include("simulation/visualization.jl")

end
