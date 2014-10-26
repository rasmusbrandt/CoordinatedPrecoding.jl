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

immutable SinglecarrierSystem <: System
    carrier_frequency::Float64
    bandwidth::Float64
end

immutable MulticarrierSystem <: System
    carrier_frequency::Float64
    bandwidth::Float64

    no_subcarriers::Int
end

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
    shadow_realization_dB::Float64
    LoS::Bool
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

    no_streams::Int

    receiver_noise_power::Float64
end

type PhysicalMS{PropagationEnvironmentState_t <: PropagationEnvironmentState} <: MS
    no_antennas::Int
    position::Position
    velocity::Velocity

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

get_no_MSs(network::Network) = length(network.MSs)
get_no_BSs(network::Network) = length(network.BSs)

get_no_antennas(node::Node) = node.no_antennas

get_no_MS_antennas(network::Network) =
    Int[ get_no_antennas(network.MSs[k]) for k = 1:get_no_MSs(network) ]
get_no_BS_antennas(network::Network) =
    Int[ get_no_antennas(network.BSs[i]) for i = 1:get_no_BSs(network) ]

get_transmit_power (BS::BS) = BS.transmit_power
set_transmit_power!(BS::BS, P::Float64) = (BS.transmit_power = P)

get_transmit_powers (network::Network) =
    Float64[ get_transmit_power(network.BSs[i]) for i = 1:get_no_BSs(network) ]
set_transmit_powers!(network::Network, P::Float64) =
    (for BS in network.BSs; set_transmit_power!(BS, P); end)

get_receiver_noise_power (MS::CanonicalMS, network::Network) =
    MS.receiver_noise_power
set_receiver_noise_power!(MS::CanonicalMS, sigma2::Float64, network::Network) =
    (MS.receiver_noise_power = sigma2)

get_receiver_noise_power (MS::PhysicalMS, network::PhysicalNetwork) =
    10^((-174 + 10*log10(network.system.bandwidth) + MS.noise_figure)/10)
set_receiver_noise_power!(MS::PhysicalMS, sigma2::Float64, network::PhysicalNetwork) =
    (MS.noise_figure = 174 + 10*log10(sigma2) - 10*log10(network.system.bandwidth))

get_receiver_noise_powers (network::Network) =
    Float64[ get_receiver_noise_power(network.MSs[k], network) for k = 1:get_no_MSs(network) ]
set_receiver_noise_powers!(network::Network, sigma2::Float64) =
    (for MS in network.MSs; set_receiver_noise_power!(MS, sigma2, network); end)

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

# ITU_R_InHNetwork
include("itu_r_inh.jl")

# Triangular3SiteNetwork
include("triangular3site.jl")
