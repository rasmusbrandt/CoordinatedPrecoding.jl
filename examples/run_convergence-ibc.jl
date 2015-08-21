#!/usr/bin/env julia

##########################################################################
# run_convergence-ic.jl
#
# Convergence as a function of number of iterations.
##########################################################################

using CoordinatedPrecoding
using Compat, JLD

##########################################################################
# Set up logging
Lumberjack.configure(Lumberjack._lumber_mill.timber_trucks["console"]; mode = "warn")
Lumberjack.add_truck(Lumberjack.LumberjackTruck("debug.log", "debug"), "debug")

##########################################################################
# General settings
srand(83196723)
start_time = Libc.strftime("%Y%m%dT%H%M%S", time())

##########################################################################
# IBC network
simulation_params = @compat Dict(
    "simulation_name" => "convergence_$(start_time)-ibc",
    "I" => 3, "Kc" => 1, "N" => 3, "M" => 3, "d" => 2,
    "Ndrops" => 10, "Nsim" => 10,
    "precoding_methods" => [
        Shi2011_WMMSE,
        Gomadam2008_MaxSINR,
        Komulainen2013_WMMSE,
        # Razaviyayn2013_MinMaxWMMSE,
        Eigenprecoding
    ],
    "aux_precoding_params" => @Compat.Dict(
        "initial_precoders" => "eigendirection",
        "stop_crit" => 0.,
        "max_iters" => 100,
    ),
    "aux_independent_variables" => [
        (set_transmit_powers_dBm!, [10, 30]),
    ],
)
network =
    setup_interfering_broadcast_channel(simulation_params["I"],
        simulation_params["Kc"], simulation_params["N"], simulation_params["M"],
        no_streams=simulation_params["d"])
raw_results = simulate_precoding_convergence(network, simulation_params)

println("-- Saving $(simulation_params["simulation_name"]) results")
save("$(simulation_params["simulation_name"]).jld",
     "simulation_params", clean_simulation_params_for_jld(simulation_params),
     "raw_results", raw_results)
