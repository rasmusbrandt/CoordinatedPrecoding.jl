#!/usr/bin/env julia

##########################################################################
# run_SNR-randomlargescale.jl
#
# Performance as a function of transmit power.
##########################################################################

using CoordinatedPrecoding
using HDF5, JLD

##########################################################################
# Set up logging
Lumberjack.configure(Lumberjack._lumber_mill.timber_trucks["console"]; mode = "warn")

##########################################################################
# General settings
srand(973472333)
start_time = strftime("%Y%m%dT%H%M%S", time())

##########################################################################
# RandomLargeScale network
simulation_params = [
    "simulation_name" => "SNR_$(start_time)-randomlargescale",
    "I" => 10, "Kc" => 1, "N" => 2, "M" => 2, "d" => 1,
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
        "stop_crit" => 1e-3,
        "max_iters" => 1000,
    ],
    "independent_variable" => (set_transmit_powers_dBm!, -50:5:0),
    # "aux_independent_variables" => [
    #     ((n, v) -> set_aux_precoding_param!(n, v, "max_iters"), [10, 50]),
    # ]
]
network =
    setup_random_large_scale_network(simulation_params["I"],
        simulation_params["Kc"], simulation_params["N"], simulation_params["M"],
        no_streams=simulation_params["d"])
raw_results, _ = simulate(network, simulation_params)

println("-- Saving $(simulation_params["simulation_name"]) results")
save("$(simulation_params["simulation_name"]).jld",
     "simulation_params", clean_simulation_params_for_jld(simulation_params),
     "raw_results", raw_results)
