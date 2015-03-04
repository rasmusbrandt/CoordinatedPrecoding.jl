#!/usr/bin/env julia

##########################################################################
# run_convergence-triangular3site.jl
#
# Convergence as a function of number of iterations.
##########################################################################

using CoordinatedPrecoding
using HDF5, JLD

##########################################################################
# Set up logging
Lumberjack.configure(Lumberjack._lumber_mill.timber_trucks["console"]; mode = "warn")
Lumberjack.add_truck(Lumberjack.LumberjackTruck("debug.log", "debug"), "debug")

##########################################################################
# General settings
srand(83196723)
start_time = strftime("%Y%m%dT%H%M%S", time())

##########################################################################
# Triangular3Site network
simulation_params = [
    "simulation_name" => "convergence_$(start_time)-triangular3site",
    "I" => 3, "Kc" => 2, "N" => 2, "M" => 4,
    "d" => 1,
    "Ndrops" => 10, "Nsim" => 10,
    "precoding_methods" => [
        Shi2011_WMMSE,
        Gomadam2008_MaxSINR,
        Komulainen2013_WMMSE,
        # Razaviyayn2013_MinMaxWMMSE,
        Eigenprecoding
    ],
    "aux_precoding_params" => [
        "initial_precoders" => "eigendirection",
        "stop_crit" => 0.,
        "max_iters" => 20,
    ],
    "aux_independent_variables" => [
        (set_transmit_powers_dBm!, [10, 30]),
    ]
]
network =
    setup_triangular3site_network(simulation_params["I"],
        simulation_params["Kc"], simulation_params["N"], simulation_params["M"],
        no_streams=simulation_params["d"])

raw_results = simulate_convergence(network, simulation_params)

println("-- Saving $(simulation_params["simulation_name"]) results")
save("$(simulation_params["simulation_name"]).jld",
     "simulation_params", clean_simulation_params_for_jld(simulation_params),
     "raw_results", raw_results)
