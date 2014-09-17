#############################################################################
# CoordinatedPrecoding
# Description here
# See http://github.com/rasmusbrandt/CoordinatedPrecoding.jl
#############################################################################

module CoordinatedPrecoding

export
# Channels
    Channel, SinglecarrierChannel, MulticarrierChannel,
    get_channel_gains_dB, get_average_channel_gains_dB,
# Networks
    Network,
    get_no_MSs, get_no_BSs, get_no_MSs_per_cell,
    get_no_antennas, get_no_MS_antennas, get_no_BS_antennas,
    get_transmit_powers, get_receiver_noise_powers, get_no_streams,
    set_transmit_powers!, set_receiver_noise_powers!, set_no_streams!,
    require_equal_no_MS_antennas, require_equal_no_BS_antennas,
    require_equal_no_streams,
    InterferenceChannel, setup_interference_channel,
    InterferingBroadcastChannel,
    Triangular3SiteNetwork, setup_triangular3site_network,
    drop_users_randomly!, draw_shadow_fading!, draw_channel,
# Cell assignment
    CellAssignment,
    get_cell_assignment, serving_BS_id, served_MS_ids,
    assign_cells_by_id, assign_cells_by_pathloss,
    assign_cells_by_instantaneous_channels, standard_cell_assignment,
# Precoding
    SinglecarrierChannel, MulticarrierChannel,
    Gomadam2008_MaxSINR, Shi2011_WMMSE, Komulainen2013_WMMSE,
    zero_receivers, initial_precoders

##########################################################################
# Channels
include("channel.jl")
# Networks
include("network/network.jl")
# Cell assignment
include("cell_assignment/cell_assignment.jl")
# Precoding algorithms
include("precoding/precoding.jl")
# Other useful stuff
include("utils.jl")
##########################################################################

end
