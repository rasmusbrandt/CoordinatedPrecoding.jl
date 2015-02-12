#!/usr/bin/env julia

##########################################################################
# run_SNR.jl
#
# Performance as a function of transmit power.
##########################################################################

using CoordinatedPrecoding
using HDF5, JLD

##########################################################################
# Set up logging
Lumberjack.configure(Lumberjack._lumber_mill.timber_trucks["console"]; mode = "warn")
Lumberjack.add_truck(Lumberjack.LumberjackTruck("debug.log", "debug"), "debug")

##########################################################################
# General settings
srand(973472333)
start_time = strftime("%Y%m%dT%H%M%S", time())

##########################################################################
# Canonical network
simulation_params = [
    "name" => "SNR_$(start_time)-ic",
    "I" => 3, "Kc" => 1, "N" => 2, "M" => 2,
    "d" => 1,
    "Ndrops" => 10, "Nsim" => 1,
    "precoding_methods" => [
        Shi2011_WMMSE,
        Gomadam2008_MaxSINR,
        Komulainen2013_WMMSE,
        Razaviyayn2013_MinMaxWMMSE,
        Eigenprecoding
    ],
    "aux_precoding_params" => [
        "initial_precoders" => "dft",
        "stop_crit" => 0.,
    ],
    "independent_variable" => (set_transmit_powers_dBm!, 0:3:30),
    "aux_independent_variables" => [
        ((n, v) -> set_aux_precoding_param!(n, v, "max_iters"), [10, 50, 100]),
    ]
]
network =
    setup_interfering_broadcast_channel(simulation_params["I"],
        simulation_params["Kc"], simulation_params["N"], simulation_params["M"],
        no_streams=simulation_params["d"])
raw_results = simulate(network, simulation_params)

println("-- Saving $(simulation_params["name"]) results")
save("$(simulation_params["name"]).jld",
     "simulation_params", clean_simulation_params_for_jld(simulation_params),
     "raw_results", raw_results)

##########################################################################
# Largescale network
simulation_params = [
    "name" => "SNR_$(start_time)-triangular3site",
    "I" => 3, "Kc" => 2, "N" => 2, "M" => 4,
    "Ps_dBm" => 0:3:30,
    "d" => 1,
    "Ndrops" => 10, "Nsim" => 1,
    "precoding_methods" => [
        Shi2011_WMMSE,
        Gomadam2008_MaxSINR,
        Komulainen2013_WMMSE,
        Razaviyayn2013_MinMaxWMMSE,
        Eigenprecoding
    ],
    "aux_precoding_params" => [
        "initial_precoders" => "dft",
        "stop_crit" => 0.,
    ],
    "independent_variable" => (set_transmit_powers_dBm!, 0:3:30),
    "aux_independent_variables" => [
        ((n, v) -> set_aux_precoding_param!(n, v, "max_iters"), [10, 50, 100]),
    ]
]
network =
    setup_triangular3site_network(simulation_params["I"],
        simulation_params["Kc"], simulation_params["N"], simulation_params["M"],
        no_streams=simulation_params["d"])
raw_results = simulate(network, simulation_params)

println("-- Saving $(simulation_params["name"]) results")
save("$(simulation_params["name"]).jld",
     "simulation_params", clean_simulation_params_for_jld(simulation_params),
     "raw_results", raw_results)
