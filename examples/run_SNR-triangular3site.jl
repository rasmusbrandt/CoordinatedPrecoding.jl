#!/usr/bin/env julia

##########################################################################
# run_SNR-triangular3site.jl
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
    "simulation_name" => "SNR_$(start_time)-triangular3site",
    "I" => 3, "Kc" => 2, "N" => 2, "M" => 4,
    "d" => 1,
    "Ndrops" => 100, "Nsim" => 1,
    "precoding_methods" => [
        Shi2011_WMMSE,
        Gomadam2008_MaxSINR,
        WeightedMaxSINR,
        Komulainen2013_WMMSE,
        # Razaviyayn2013_MinMaxWMMSE,
        Eigenprecoding
    ],
    "aux_precoding_params" => [
        "initial_precoders" => "eigendirection",
        "stop_crit" => 1e-3,
        "max_iters" => 1000,
    ],
    "independent_variable" => (set_transmit_powers_dBm!, -10:3:30),
    # "aux_independent_variables" => [
    #     ((n, v) -> set_aux_precoding_param!(n, v, "max_iters"), [10, 50]),
    # ]
]
network =
    setup_triangular3site_network(simulation_params["I"],
        simulation_params["Kc"], simulation_params["N"], simulation_params["M"],
        no_streams=simulation_params["d"])
raw_results = simulate(network, simulation_params)

println("-- Saving $(simulation_params["simulation_name"]) results")
save("$(simulation_params["simulation_name"]).jld",
     "simulation_params", clean_simulation_params_for_jld(simulation_params),
     "raw_results", raw_results)
