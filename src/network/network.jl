##########################################################################
# Geography
immutable Position
    x::Float64
    y::Float64
end

+(A::Position, B::Position) = Position(A.x + B.x, A.y + B.y)
-(A::Position, B::Position) = Position(A.x - B.x, A.y - B.y)

get_distance(A::Position, B::Position) = sqrt((A.x - B.x)^2 + (A.y - B.y)^2)
get_angle(A::Position, B::Position) = atan2(A.y - B.y, A.x - B.x)

immutable Velocity
    x::Float64
    y::Float64
end

##########################################################################
# Systems
abstract System

typealias AuxPrecodingParams Dict{ASCIIString, Any} # belongs in precoding.jl
typealias AuxCellAssignmentParams Dict{ASCIIString, Any} # belongs in cell_assignment.jl
typealias AuxClusterAssignmentParams Dict{ASCIIString, Any} # belongs in cluster_assignment.jl

type SinglecarrierSystem <: System
    carrier_frequency::Float64
    bandwidth::Float64

    aux_precoding_params::AuxPrecodingParams
    aux_cell_assignment_params::AuxCellAssignmentParams
    aux_cluster_assignment_params::AuxClusterAssignmentParams
end
SinglecarrierSystem() =
    SinglecarrierSystem(0, 0, AuxPrecodingParams(), AuxCellAssignmentParams(), AuxClusterAssignmentParams())
SinglecarrierSystem(carrier_frequency::Float64, bandwidth::Float64) =
    SinglecarrierSystem(carrier_frequency, bandwidth, AuxPrecodingParams(), AuxCellAssignmentParams(), AuxClusterAssignmentParams())

type MulticarrierSystem <: System
    carrier_frequency::Float64
    bandwidth::Float64

    no_subcarriers::Int

    aux_precoding_params::AuxPrecodingParams
    aux_cell_assignment_params::AuxCellAssignmentParams
    aux_cluster_assignment_params::AuxClusterAssignmentParams
end
MulticarrierSystem() =
    MulticarrierSystem(0, 0, 1, AuxPrecodingParams(), AuxCellAssignmentParams(), AuxClusterAssignmentParams())
MulticarrierSystem(carrier_frequency::Float64, bandwidth::Float64, no_subcarriers::Int) =
    MulticarrierSystem(carrier_frequency, bandwidth, no_subcarriers, AuxPrecodingParams(), AuxCellAssignmentParams(), AuxClusterAssignmentParams())

##########################################################################
# Antenna params
abstract AntennaParams

immutable OmnidirectionalAntennaParams <: AntennaParams
    antenna_gain_dB::Float64
end

get_antenna_gain(antenna_params::OmnidirectionalAntennaParams) = 10^(antenna_params.antenna_gain_dB/10)
get_antenna_gain(antenna_params::OmnidirectionalAntennaParams, angle::Float64) = get_antenna_gain(antenna_params)

##########################################################################
# Propagation environments
abstract PropagationEnvironment
abstract PropagationEnvironmentState

immutable SimpleLargescaleFadingEnvironment <: PropagationEnvironment
    pathloss_alpha::Float64
    pathloss_beta::Float64
    penetration_loss_dB::Float64
    shadow_sigma_dB::Float64
end

immutable SimpleLargescaleFadingEnvironmentState <: PropagationEnvironmentState
    shadow_realizations_dB::Vector{Float64}
    LoSs::BitArray
end

##########################################################################
# Nodes
abstract Node

abstract BS <: Node

type CanonicalBS <: BS
    no_antennas::Int

    transmit_power::Float64
end

type PhysicalBS{AntennaParams_t <: AntennaParams} <: BS
    no_antennas::Int
    position::Position

    transmit_power::Float64

    antenna_params::AntennaParams_t
end

abstract MS <: Node

type CanonicalMS <: MS
    no_antennas::Int

    user_priority::Float64
    no_streams::Int

    receiver_noise_power::Float64
end

type PhysicalMS{PropagationEnvironmentState_t <: PropagationEnvironmentState} <: MS
    no_antennas::Int
    position::Position
    velocity::Velocity

    user_priority::Float64
    no_streams::Int

    antenna_gain_dB::Float64
    noise_figure::Float64
    propagation_environment_state::PropagationEnvironmentState_t
end

##########################################################################
# Networks
abstract Network
abstract CanonicalNetwork <: Network
abstract PhysicalNetwork <: Network

get_cell_assignment (network::Network) = network.cell_assignment
set_cell_assignment!(network::Network, cell_assignment::CellAssignment) =
    (network.cell_assignment = cell_assignment)

get_cluster_assignment (network::Network) = network.cluster_assignment
set_cluster_assignment!(network::Network, cluster_assignment::ClusterAssignment) =
    (network.cluster_assignment = cluster_assignment)

# FIXME: Clean this up using macros
get_aux_precoding_param (network::Network, k::ASCIIString) =
    (network.system.aux_precoding_params[k])
set_aux_precoding_param!(network::Network, v, k::ASCIIString) =
    (network.system.aux_precoding_params[k] = v)
get_aux_precoding_params (network::Network) = network.system.aux_precoding_params
set_aux_precoding_params!(network::Network, additional::AuxPrecodingParams) =
    merge!(network.system.aux_precoding_params, additional)

get_aux_cell_assignment_param (network::Network, k::ASCIIString) =
    (network.system.aux_cell_assignment_params[k])
set_aux_cell_assignment_param!(network::Network, v, k::ASCIIString) =
    (network.system.aux_cell_assignment_params[k] = v)
get_aux_cell_assignment_params (network::Network) = network.system.aux_cell_assignment_params
set_aux_cell_assignment_params!(network::Network, additional::AuxCellAssignmentParams) =
    merge!(network.system.aux_cell_assignment_params, additional)

get_aux_cluster_assignment_param (network::Network, k::ASCIIString) =
    (network.system.aux_cluster_assignment_params[k])
set_aux_cluster_assignment_param!(network::Network, v, k::ASCIIString) =
    (network.system.aux_cluster_assignment_params[k] = v)
get_aux_cluster_assignment_params (network::Network) = network.system.aux_cluster_assignment_params
set_aux_cluster_assignment_params!(network::Network, additional::AuxClusterAssignmentParams) =
    merge!(network.system.aux_cluster_assignment_params, additional)

get_no_MSs(network::Network) = length(network.MSs)
get_no_BSs(network::Network) = length(network.BSs)

get_no_antennas(node::Node) = node.no_antennas

get_no_MS_antennas(network::Network) =
    Int[ get_no_antennas(network.MSs[k]) for k = 1:get_no_MSs(network) ]
get_no_BS_antennas(network::Network) =
    Int[ get_no_antennas(network.BSs[i]) for i = 1:get_no_BSs(network) ]

get_transmit_power (BS::BS) = BS.transmit_power
set_transmit_power!(BS::BS, P::Float64) = (BS.transmit_power = P)
get_transmit_power_dBm (BS::BS) = 10*log10(get_transmit_power(BS))
set_transmit_power_dBm!(BS::BS, PdBm) = set_transmit_power!(BS, 10^(PdBm/10))

get_transmit_powers (network::Network) =
    Float64[ get_transmit_power(network.BSs[i]) for i = 1:get_no_BSs(network) ]
set_transmit_powers!(network::Network, P) =
    (for BS in network.BSs; set_transmit_power!(BS, P); end)
get_transmit_powers_dBm (network::Network) =
    Float64[ get_transmit_power_dBm(network.BSs[i]) for i = 1:get_no_BSs(network) ]
set_transmit_powers_dBm!(network::Network, PdBm) =
    (for BS in network.BSs; set_transmit_power_dBm!(BS, PdBm); end)

get_receiver_noise_power (MS::CanonicalMS, network::Network) =
    MS.receiver_noise_power
set_receiver_noise_power!(MS::CanonicalMS, sigma2, network::Network) =
    (MS.receiver_noise_power = sigma2)
get_receiver_noise_power_dBm (MS::CanonicalMS, network::Network) =
    10*log10(get_receiver_noise_power(MS, network))
set_receiver_noise_power_dBm!(MS::CanonicalMS, sigma2dBm, network::Network) =
    set_receiver_noise_power(MS, 10^(sigma2dBm/10), network)

get_receiver_noise_power (MS::PhysicalMS, network::PhysicalNetwork) =
    10^(get_receiver_noise_power_dBm(MS, network)/10)
set_receiver_noise_power!(MS::PhysicalMS, sigma2, network::PhysicalNetwork) =
    (MS.noise_figure = 174. + 10*log10(sigma2) - 10*log10(network.system.bandwidth))
get_receiver_noise_power_dBm (MS::PhysicalMS, network::PhysicalNetwork) =
    (-174. + 10*log10(network.system.bandwidth) + MS.noise_figure)
set_receiver_noise_power_dBm!(MS::PhysicalMS, sigma2dBm, network::PhysicalNetwork) =
    (MS.noise_figure = 174. + sigma2dBm - 10*log10(network.system.bandwidth))

get_receiver_noise_powers (network::Network) =
    Float64[ get_receiver_noise_power(network.MSs[k], network) for k = 1:get_no_MSs(network) ]
set_receiver_noise_powers!(network::Network, sigma2) =
    (for MS in network.MSs; set_receiver_noise_power!(MS, sigma2, network); end)
get_receiver_noise_powers_dBm (network::Network) =
    Float64[ get_receiver_noise_power_dBm(network.MSs[k], network) for k = 1:get_no_MSs(network) ]
set_receiver_noise_powers_dBm!(network::Network, sigma2dBm) =
    (for MS in network.MSs; set_receiver_noise_power_dBm!(MS, sigma2dBm, network); end)

get_user_priority (MS::MS) = MS.user_priority
set_user_priority!(MS::MS, α) = (MS.user_priority = α)

get_user_priorities (network::Network) =
    Float64[ get_user_priority(network.MSs[k]) for k = 1:get_no_MSs(network) ]
set_user_priorities!(network::Network, α) =
    (for MS in network.MSs; set_user_priority(MS.no_streams, α); end)

get_no_streams (MS::MS) = MS.no_streams
set_no_streams!(MS::MS, d::Int) = (MS.no_streams = d)

get_no_streams (network::Network) =
    Int[ get_no_streams(network.MSs[k]) for k = 1:get_no_MSs(network) ]
set_no_streams!(network::Network, d::Int) =
    (for MS in network.MSs; set_no_streams(MS.no_streams, d); end)

require_equal_no_MS_antennas(network::Network) =
    (Ns = get_no_MS_antennas(network); all(Ns .== Ns[1]) || error("MSs must all have the same number of antennas."))
require_equal_no_BS_antennas(network::Network) =
    (Ms = get_no_BS_antennas(network); all(Ms .== Ms[1]) || error("BSs must all have the same number of antennas."))

require_single_antenna_MSs(network::Network) =
    (Ns = get_no_MS_antennas(network); all(Ns .== 1) || error("MSs must not have multiple antennas."))
require_single_antenna_BSs(network::Network) =
    (Ms = get_no_BS_antennas(network); all(Ms .== 1) || error("BSs must not have multiple antennas."))

require_equal_no_streams(network::Network) =
    (ds = get_no_streams(network); all(ds .== ds[1]) || error("MSs must all have the same number of streams."))

require_single_stream(network::Network) =
    (ds = get_no_streams(network); all(ds .== 1) || error("MSs must not have multiple streams."))

get_distances(network::PhysicalNetwork) = 
    Float64[ get_distance(network.MSs[k].position, network.BSs[i].position) for k = 1:get_no_MSs(network), i = 1:get_no_BSs(network) ]

get_angles(network::PhysicalNetwork) = 
    Float64[ (get_angle(network.MSs[k].position, network.BSs[i].position) - network.BSs[i].antenna_params.bore_sight_angle) for k = 1:get_no_MSs(network), i = 1:get_no_BSs(network) ]

# Interfering broadcast channel
include("ibc.jl")

# IndoorsNetwork
include("indoors.jl")

# Triangular3SiteNetwork
include("triangular3site.jl")
