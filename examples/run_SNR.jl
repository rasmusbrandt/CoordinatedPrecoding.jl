#!/usr/bin/env julia

##########################################################################
# run_SNR.jl
#
# Performance as a function of transmit power.
##########################################################################

using CoordinatedPrecoding
using HDF5, JLD

##########################################################################
# General settings
srand(973472333)
start_time = strftime("%Y%m%dT%H%M%S", time())

##########################################################################
# Canonical network
simulation_params = [
    "name" => "$(start_time)-canonical",
    "K" => 3, "N" => 2, "M" => 2,
    "Ps_dBm" => 0:3:30,
    "d" => 1,
    "Ndrops" => 10, "Nsim" => 5,
    "precoding_methods" => [
        Shi2011_WMMSE,
        Gomadam2008_MaxSINR,
        Komulainen2013_WMMSE,
        Razaviyayn2013_MaxMinWMMSE,
        Eigenprecoding
    ]
]
precoding_settings = [
    "stop_crit" => 20,
    "initial_precoders" => "dft",
]
network =
    setup_interference_channel(simulation_params["K"],
        simulation_params["N"], simulation_params["M"],
        no_streams=simulation_params["d"])

results = simulate_SNR(network, simulation_params, precoding_settings)

println("-- Saving $(simulation_params["name"]) results")
save("SNR_$(simulation_params["name"]).jld",
     "simulation_params", clean_simulation_params_for_JLD(simulation_params),
     "precoding_settings", precoding_settings,
     "results", results)

##########################################################################
# Largescale network
simulation_params = [
    "name" => "$(start_time)-largescale",
    "Kc" => 2, "N" => 2, "M" => 4,
    "Ps_dBm" => 0:3:30,
    "d" => 1,
    "Ndrops" => 1, "Nsim" => 1,
    "precoding_methods" => [
        Shi2011_WMMSE,
        Gomadam2008_MaxSINR,
        Komulainen2013_WMMSE,
        Razaviyayn2013_MaxMinWMMSE,
        Eigenprecoding
    ]
]
precoding_settings = [
    "stop_crit" => 20,
    "initial_precoders" => "dft",
]
network =
    setup_triangular3site_network(simulation_params["Kc"],
        simulation_params["N"], simulation_params["M"],
        no_streams=simulation_params["d"])

results = simulate_SNR(network, simulation_params, precoding_settings)

println("-- Saving $(simulation_params["name"]) results")
save("SNR_$(simulation_params["name"]).jld",
     "simulation_params", clean_simulation_params_for_JLD(simulation_params),
     "precoding_settings", precoding_settings,
     "results", results)