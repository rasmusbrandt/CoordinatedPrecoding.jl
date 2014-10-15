#!/usr/bin/env julia

##########################################################################
# plot_convergence.jl
#
# Plots convergence curves.
##########################################################################

##########################################################################
# Load data
#
# Do this before loading other code, otherwise the JLD module might crash!
using HDF5, JLD, ArgParse
s = ArgParseSettings()
@add_arg_table s begin
    "file_name"
        help = "file name with results"
        required = true
end
parsed_args = parse_args(s)
data = load(parsed_args["file_name"])

##########################################################################
# Plot parameters
using CoordinatedPrecoding

plot_params = [
    "figsize" => (8,4),
    "precoding_methods" =>
        [
          "Shi2011_WMMSE" => [ ("rates", "b-", "WMMSE bound"), ],
          "Gomadam2008_MaxSINR" => [ ("rates", "r-", "MaxSINR bound"), ],
          "Razaviyayn2013_MaxMinWMMSE" => [ ("rates", "b-.", "MaxMinWMMSE bound"), ],
          "Komulainen2013_WMMSE" => [ ("rates", "m-", "Rounded bound"), ],
          "Eigenprecoding" =>
            [ ("intercell_tdma_rates", "c-", "TDMA bound"),
              ("intracell_tdma_rates", "c-.", "Intracell TDMA bound"),
              ("uncoord_rates", "k-", "Uncoordinated bound"),
            ],
        ]
]

##########################################################################
# Plot it
plot_convergence(
    data["results"],
    data["simulation_params"],
    data["precoding_settings"],
    plot_params)
