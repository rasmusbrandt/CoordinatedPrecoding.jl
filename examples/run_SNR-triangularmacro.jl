#!/usr/bin/env julia

##########################################################################
# run_SNR-triangularmacro.jl
#
# Performance as a function of transmit power.
##########################################################################

using CoordinatedPrecoding
using Compat, JLD

##########################################################################
# Set up logging
Lumberjack.configure(Lumberjack._lumber_mill.timber_trucks["console"]; mode = "warn")

##########################################################################
# General settings
srand(973472333)
start_time = Libc.strftime("%Y%m%dT%H%M%S", time())

##########################################################################
# TriangularMacro network
simulation_params = @compat Dict(
    "simulation_name" => "SNR_$(start_time)-triangularmacro",
    "Kc" => 2, "N" => 2, "M" => 4, "d" => 1,
    "Ndrops" => 10, "Nsim" => 10,
    "precoding_methods" => [
        Shi2011_WMMSE,
        Gomadam2008_MaxSINR,
        Komulainen2013_WMMSE,
        # Razaviyayn2013_MinMaxWMMSE,
        Eigenprecoding
    ],
    "aux_precoding_params" => Dict(
        "initial_precoders" => "eigendirection",
        "stop_crit" => 1e-3,
        "max_iters" => 1000,
    ),
    "independent_variable" => (set_transmit_powers_dBm!, -20:5:30),
    # "aux_independent_variables" => [
    #     ((n, v) -> set_aux_precoding_param!(n, v, "max_iters"), [10, 50]),
    # ],
)
network =
    setup_triangularmacro_network(simulation_params["Kc"],
        simulation_params["N"], simulation_params["M"],
        num_streams=simulation_params["d"])
raw_results, _ = simulate(network, simulation_params)

println("-- Saving $(simulation_params["simulation_name"]) results")
save("$(simulation_params["simulation_name"]).jld",
     "simulation_params", clean_simulation_params_for_jld(simulation_params),
     "raw_results", raw_results)
