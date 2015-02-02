#!/usr/bin/env julia

##########################################################################
# plot_SNR.jl
#
# Plots SNR curves.
##########################################################################

using CoordinatedPrecoding
using HDF5, JLD, ArgParse

##########################################################################
# Load data
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
plot_params = [
    "name_suffix" => "",

    "figsize" => (8,5),

    "objectives" => [
        "sumrate" => (r -> sum(r, 5:6), [ "xlabel" => "Transmit power [dBm]", "ylabel" => "Sum rate [bits/s/Hz]" ]),
        "minrate" => (r -> minimum(sum(r, 6), 5), [ "xlabel" => "Transmit power [dBm]", "ylabel" => "Min rate [bits/s/Hz]", ]),
    ],

    "precoding_methods" => {
        "Shi2011_WMMSE" => [
            ("logdet_rates", [ "key" => "b-", "legend" => "WMMSE" ]),
            ("MMSE_rates", [ "key" => "b--", "legend" => "WMMSE" ]),
        ],

        "Gomadam2008_MaxSINR" => [
            ("logdet_rates", [ "key" => "r-", "legend" => "MaxSINR" ]),
            ("MMSE_rates", [ "key" => "r--", "legend" => "MaxSINR" ]),
        ],

        "Razaviyayn2013_MinMaxWMMSE" => [
            ("logdet_rates", [ "key" => "g-", "legend" => "MinMax-WMMSE" ]),
            ("MMSE_rates", [ "key" => "g--", "legend" => "MinMax-WMMSE" ]),
        ],

        "Komulainen2013_WMMSE" => [
            ("logdet_rates", [ "key" => "m-", "legend" => "DiagWMMSE" ]),
            ("MMSE_rates", [ "key" => "m--", "legend" => "DiagWMMSE" ]),
        ],

        "Eigenprecoding" => {
            ("intercell_tdma_logdet_rates", [ "key" => "c-", "legend" => "TDMA" ]),
            ("intracell_tdma_logdet_rates", [ "key" => "c-.", "legend" => "Intracell TDMA" ]),
            ("uncoord_logdet_rates", [ "key" => "k-", "legend" => "Uncoordinated" ]),
        },
    },
]

##########################################################################
# Plot it
processed_results = process(data["raw_results"], data["simulation_params"], plot_params)
plot(processed_results, data["simulation_params"], plot_params)
