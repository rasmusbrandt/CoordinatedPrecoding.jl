##########################################################################
# Geography
immutable Position
    x::Float64
    y::Float64
end

+(A::Position, B::Position) = Position(A.x + B.x, A.y + B.y)
-(A::Position, B::Position) = Position(A.x - B.x, A.y - B.y)
rotate(p::Position, θ) = Position(cos(θ)*p.x - sin(θ)*p.y, sin(θ)*p.x + cos(θ)*p.y)

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
typealias AuxAssignmentParams Dict{ASCIIString, Any} # belongs in assignment.jl

type SinglecarrierSystem <: System
    carrier_frequency::Float64
    bandwidth::Float64

    aux_precoding_params::AuxPrecodingParams
    aux_assignment_params::AuxAssignmentParams
end
SinglecarrierSystem() =
    SinglecarrierSystem(0, 0, AuxPrecodingParams(), AuxAssignmentParams())
SinglecarrierSystem(carrier_frequency, bandwidth) =
    SinglecarrierSystem(carrier_frequency, bandwidth, AuxPrecodingParams(), AuxAssignmentParams())

type MulticarrierSystem <: System
    carrier_frequency::Float64
    bandwidth::Float64

    no_subcarriers::Int

    aux_precoding_params::AuxPrecodingParams
    aux_assignment_params::AuxAssignmentParams
end
MulticarrierSystem() =
    MulticarrierSystem(0, 0, 1, AuxPrecodingParams(), AuxAssignmentParams())
MulticarrierSystem(carrier_frequency, bandwidth, no_subcarriers) =
    MulticarrierSystem(carrier_frequency, bandwidth, no_subcarriers, AuxPrecodingParams(), AuxAssignmentParams())

##########################################################################
# Antenna params
abstract AntennaParams

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

get_no_antennas(node) = node.no_antennas

get_user_priority (MS) = MS.user_priority
set_user_priority!(MS, α) = (MS.user_priority = α)

get_no_streams (MS::MS) = MS.no_streams
set_no_streams!(MS::MS, d) = (MS.no_streams = d)

get_receiver_noise_power (MS::CanonicalMS, network) = MS.receiver_noise_power
set_receiver_noise_power!(MS::CanonicalMS, sigma2, network) = (MS.receiver_noise_power = sigma2)
get_receiver_noise_power_dBm (MS::CanonicalMS, network) = 10*log10(get_receiver_noise_power(MS, network))
set_receiver_noise_power_dBm!(MS::CanonicalMS, sigma2_dBm, network) = set_receiver_noise_power(MS, 10^(sigma2_dBm/10), network)

get_receiver_noise_power (MS::PhysicalMS, network) = 10^(get_receiver_noise_power_dBm(MS, network)/10)
set_receiver_noise_power!(MS::PhysicalMS, sigma2, network) = (MS.noise_figure = 174. + 10*log10(sigma2) - 10*log10(network.system.bandwidth))
get_receiver_noise_power_dBm (MS::PhysicalMS, network) = (-174. + 10*log10(network.system.bandwidth) + MS.noise_figure)
set_receiver_noise_power_dBm!(MS::PhysicalMS, sigma2_dBm, network) = (MS.noise_figure = 174. + sigma2_dBm - 10*log10(network.system.bandwidth))

get_transmit_power (BS) = BS.transmit_power
set_transmit_power!(BS, P) = (BS.transmit_power = P)
get_transmit_power_dBm (BS) = 10*log10(get_transmit_power(BS))
set_transmit_power_dBm!(BS, P_dBm) = set_transmit_power!(BS, 10^(P_dBm/10))

get_distance(MS, BS) = get_distance(MS.position, BS.position)

##########################################################################
# OmnidirectionalAntennaParams (not angle dependent)
immutable OmnidirectionalAntennaParams <: AntennaParams
    antenna_gain_dB::Float64
end

get_antenna_gain(antenna_params::OmnidirectionalAntennaParams) =
    10^(antenna_params.antenna_gain_dB/10)
get_antenna_gain(antenna_params::OmnidirectionalAntennaParams, angle) =
    get_antenna_gain(antenna_params)

get_angle{AntennaParams_t <: OmnidirectionalAntennaParams}(MS, BS::PhysicalBS{AntennaParams_t}) =
    get_angle(MS.position, BS.position)

##########################################################################
# Networks
abstract Network
abstract CanonicalNetwork <: Network
abstract PhysicalNetwork <: Network

get_no_MSs(network) = length(network.MSs)
get_no_BSs(network) = length(network.BSs)

get_assignment (network) = network.assignment
set_assignment!(network, assignment) = (network.assignment = assignment)

get_aux_precoding_param (network, k) = (network.system.aux_precoding_params[k])
set_aux_precoding_param!(network, v, k) = (network.system.aux_precoding_params[k] = v)
get_aux_precoding_params (network) = network.system.aux_precoding_params
set_aux_precoding_params!(network, additional) = merge!(network.system.aux_precoding_params, additional)

get_aux_assignment_param (network, k) = (network.system.aux_assignment_params[k])
set_aux_assignment_param!(network, v, k) = (network.system.aux_assignment_params[k] = v)
get_aux_assignment_params (network) = network.system.aux_assignment_params
set_aux_assignment_params!(network, additional) = merge!(network.system.aux_assignment_params, additional)

get_no_MS_antennas(network) = [ get_no_antennas(network.MSs[k]) for k = 1:get_no_MSs(network) ]
get_no_BS_antennas(network) = [ get_no_antennas(network.BSs[i]) for i = 1:get_no_BSs(network) ]

get_transmit_powers (network) = [ get_transmit_power(network.BSs[i]) for i = 1:get_no_BSs(network) ]
set_transmit_powers!(network, P::Real) = (for BS in network.BSs; set_transmit_power!(BS, P); end)
set_transmit_powers!(network, Ps::Vector) = (for i = 1:length(Ps); set_transmit_power!(network.BSs[i], Ps[i]); end)
get_transmit_powers_dBm (network) = [ get_transmit_power_dBm(network.BSs[i]) for i = 1:get_no_BSs(network) ]
set_transmit_powers_dBm!(network, P_dBm::Real) = (for BS in network.BSs; set_transmit_power_dBm!(BS, P_dBm); end)
set_transmit_powers_dBm!(network, Ps_dBm::Vector) = (for i = 1:length(Ps_dBm); set_transmit_power_dBm!(network.BSs[i], Ps_dBm[i]); end)

get_receiver_noise_powers (network) = [ get_receiver_noise_power(network.MSs[k], network) for k = 1:get_no_MSs(network) ]
set_receiver_noise_powers!(network, sigma2::Real) = (for MS in network.MSs; set_receiver_noise_power!(MS, sigma2, network); end)
set_receiver_noise_powers!(network, sigma2s::Vector) = (for k = 1:length(sigma2s); set_receiver_noise_power!(network.MSs[k], sigma2s[k], network); end)
get_receiver_noise_powers_dBm (network) = [ get_receiver_noise_power_dBm(network.MSs[k], network) for k = 1:get_no_MSs(network) ]
set_receiver_noise_powers_dBm!(network, sigma2_dBm::Real) = (for MS in network.MSs; set_receiver_noise_power_dBm!(MS, sigma2_dBm, network); end)
set_receiver_noise_powers_dBm!(network, sigma2s_dBm::Vector) = (for k = 1:length(sigma2s_dBm); set_receiver_noise_power_dBm!(network.MSs[k], sigma2s_dBm[k], network); end)

get_user_priorities (network) = [ get_user_priority(network.MSs[k]) for k = 1:get_no_MSs(network) ]
set_user_priorities!(network, α::Real) = (for MS in network.MSs; set_user_priority(MS, α); end)
set_user_priorities!(network, α::Vector) = (for k = 1:length(α); set_user_priority(network.MSs[k], α[k]); end)

get_no_streams (network::Network) = [ get_no_streams(network.MSs[k]) for k = 1:get_no_MSs(network) ]
set_no_streams!(network::Network, d::Int) = (for MS in network.MSs; set_no_streams(MS, d); end)
set_no_streams!(network::Network, ds::Vector) = (for k = 1:length(ds); set_no_streams(network.MSs[k], ds[k]); end)

require_equal_no_MS_antennas(network) = (Ns = get_no_MS_antennas(network); all(Ns .== Ns[1]) || Lumberjack.error("MSs must all have the same number of antennas."))
require_equal_no_BS_antennas(network) = (Ms = get_no_BS_antennas(network); all(Ms .== Ms[1]) || Lumberjack.error("BSs must all have the same number of antennas."))
require_single_antenna_MSs(network) = (Ns = get_no_MS_antennas(network); all(Ns .== 1) || Lumberjack.error("MSs must not have multiple antennas."))
require_single_antenna_BSs(network) = (Ms = get_no_BS_antennas(network); all(Ms .== 1) || Lumberjack.error("BSs must not have multiple antennas."))
require_equal_no_streams(network) = (ds = get_no_streams(network); all(ds .== ds[1]) || Lumberjack.error("MSs must all have the same number of streams."))
require_single_stream(network) = (ds = get_no_streams(network); all(ds .== 1) || Lumberjack.error("MSs must not have multiple streams."))

get_distances(network::PhysicalNetwork) = [ get_distance(network.MSs[k], network.BSs[i]) for k = 1:get_no_MSs(network), i = 1:get_no_BSs(network) ]
get_angles(network::PhysicalNetwork) = [ get_angle(network.MSs[k], network.BSs[i]) for k = 1:get_no_MSs(network), i = 1:get_no_BSs(network) ]

##########################################################################
# Include network implementations
include("IntereringBroadcastChannel.jl")
include("IndoorsNetwork.jl")
include("Triangular3SiteNetwork.jl")
include("TriangularHetNetNetwork.jl")
