#!/usr/bin/env julia

##########################################################################
# run_performancetest.jl
#
# Compares runtimes of different methods.
##########################################################################

using CoordinatedPrecoding
using HDF5, JLD

##########################################################################
# General settings
srand(8071232234)
start_time = strftime("%Y%m%dT%H%M%S", time())

##########################################################################
# Performance test
simulation_params = {
    "K" => 3, "N" => 2, "M" => 2,
    "P_dBm" => 20.,
    "d" => 1,
    "Ntest" => 10,
    "precoding_methods" => [
        Shi2011_WMMSE,
        Gomadam2008_MaxSINR,
        Komulainen2013_WMMSE,
        Razaviyayn2013_MaxMinWMMSE,
        Eigenprecoding
    ]
}
precoding_settings = {
    "stop_crit" => 20,
    "initial_precoders" => "dft",
}
network =
    setup_interference_channel(simulation_params["K"],
        simulation_params["N"], simulation_params["M"],
        transmit_power=10^(simulation_params["P_dBm"]/10),
        no_streams=simulation_params["d"])

perform_performancetest(network, simulation_params, precoding_settings)
