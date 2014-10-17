#!/usr/bin/env julia

##########################################################################
# run_convergence.jl
#
# Convergence as a function of number of iterations.
##########################################################################

using CoordinatedPrecoding
using HDF5, JLD

##########################################################################
# General settings
srand(83196723)
start_time = strftime("%Y%m%dT%H%M%S", time())

##########################################################################
# Canonical network
simulation_params = [
    "name" => "$(start_time)-canonical",
    "I" => 3, "Kc" => 1, "N" => 2, "M" => 2,
    "P_dBm" => 20.,
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
    "stop_crit" => 100,
    "initial_precoders" => "dft",
]
network =
    setup_interfering_broadcast_channel(simulation_params["I"],
        simulation_params["Kc"], simulation_params["N"], simulation_params["M"],
        transmit_power=10^(simulation_params["P_dBm"]/10),
        no_streams=simulation_params["d"])

results = simulate_convergence(network, simulation_params, precoding_settings)

println("-- Saving $(simulation_params["name"]) results")
save("convergence_$(simulation_params["name"]).jld",
     "simulation_params", clean_simulation_params_for_JLD(simulation_params),
     "precoding_settings", precoding_settings,
     "results", results)

##########################################################################
# Largescale network
simulation_params = [
    "name" => "$(start_time)-largescale",
    "I" => 3, "Kc" => 2, "N" => 2, "M" => 4,
    "P_dBm" => 18.2,
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
    "stop_crit" => 100,
    "initial_precoders" => "dft",
]
network =
    setup_triangular3site_network(simulation_params["I"],
        simulation_params["Kc"], simulation_params["N"], simulation_params["M"],
        transmit_power=10^(simulation_params["P_dBm"]/10),
        no_streams=simulation_params["d"])

results = simulate_convergence(network, simulation_params, precoding_settings)

println("-- Saving $(simulation_params["name"]) results")
save("convergence_$(simulation_params["name"]).jld",
     "simulation_params", clean_simulation_params_for_JLD(simulation_params),
     "precoding_settings", precoding_settings,
     "results", results)
