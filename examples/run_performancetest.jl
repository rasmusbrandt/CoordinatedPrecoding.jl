#!/usr/bin/env julia

##########################################################################
# run_performancetest.jl
#
# Compares runtimes of different methods.
##########################################################################

using CoordinatedPrecoding

##########################################################################
# Set up logging
Lumberjack.configure(Lumberjack._lumber_mill.timber_trucks["console"]; mode = "warn")
Lumberjack.add_truck(Lumberjack.LumberjackTruck("debug.log", "debug"), "debug")

##########################################################################
# General settings
srand(8071232234)
start_time = strftime("%Y%m%dT%H%M%S", time())

##########################################################################
# Performance test
simulation_params = [
    "I" => 3, "Kc" => 1, "N" => 2, "M" => 2,
    "P_dBm" => 20.,
    "d" => 1,
    "Ntest" => 100,
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
        "max_iters" => 100,
    ]
]
network =
    setup_interfering_broadcast_channel(simulation_params["I"],
        simulation_params["Kc"], simulation_params["N"], simulation_params["M"],
        transmit_power=10^(simulation_params["P_dBm"]/10),
        no_streams=simulation_params["d"])
set_aux_precoding_params!(network, simulation_params["aux_precoding_params"])

simulate_performance(network, simulation_params)
