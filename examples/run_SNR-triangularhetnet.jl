#!/usr/bin/env julia

##########################################################################
# run_SNR-triangularhetnet.jl
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
# Triangular3Site network
simulation_params = [
    "simulation_name" => "SNR_$(start_time)-triangularhetnet",
    "I" => 1, "Kc" => 2, "N" => 2, "M" => 2,
    "d" => 1,
    "Ndrops" => 10, "Nsim" => 10,
    "assignment_methods" => [
        assign_cells_by_large_scale_fading!,
    ],
    "precoding_methods" => [
        Shi2011_WMMSE,
        Gomadam2008_MaxSINR,
        Komulainen2013_WMMSE,
        # Razaviyayn2013_MinMaxWMMSE,
        Eigenprecoding
    ],
    "aux_precoding_params" => [
        "initial_precoders" => "eigendirection",
        "stop_crit" => 1e-3,
        "max_iters" => 1000,
    ],
    "independent_variable" => (set_transmit_powers_dBm!, [ p*ones(6) - [0,0,0,20,20,20] for p = -20:5:30 ]),
    # "aux_independent_variables" => [
    #     ((n, v) -> set_aux_precoding_param!(n, v, "max_iters"), [10, 50]),
    # ]
]
network =
    setup_triangularhetnet_network(simulation_params["I"],
        simulation_params["Kc"], simulation_params["N"], simulation_params["M"],
        no_streams=simulation_params["d"])
raw_results = simulate(network, simulation_params)

println("-- Saving $(simulation_params["simulation_name"]) results")
save("$(simulation_params["simulation_name"]).jld",
     "simulation_params", clean_simulation_params_for_jld(simulation_params),
     "raw_results", raw_results)
