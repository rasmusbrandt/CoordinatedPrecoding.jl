#!/usr/bin/env julia

##########################################################################
# plot_convergence.jl
#
# Plots convergence curves.
##########################################################################

using CoordinatedPrecoding
using HDF5, JLD, ArgParse

##########################################################################
# Load data
s = ArgParseSettings()
@add_arg_table s begin
    "file_names"
        help = "file name with results"
        required = true
        nargs = '+'
end
parsed_args = parse_args(s)

##########################################################################
# Plot parameters
plot_params = [
    "plot_name" => "",

    "objective" => :sumrate,

    "figure" => [
        :figsize => (8,5),
        :dpi => 125,
        :facecolor => "w",
        :edgecolor => "k",
    ],

    "axes" => [
        :xlabel => "Transmit power [dBm]",
        :ylabel => "Sum rate [bits/s/Hz]",
        :xscale => "linear",
        :yscale => "linear",
    ],

    "legend" => [
        :loc => "best",
        :fontsize => 8,
    ],

    # "confidence_interval_zalpha2" => 1.96,

    "precoding_methods" => [
        "Shi2011_WMMSE" => [
            ("logdet_rates", [ :color => "b", :linestyle => "-", :label => "WMMSE (logdet)" ]),
            ("MMSE_rates", [ :color => "b", :linestyle => "--",  :label => "WMMSE (MMSE)" ]),
        ],

        "Gomadam2008_MaxSINR" => [
            ("logdet_rates", [ :color => "r", :linestyle => "-", :label => "MaxSINR (logdet)" ]),
            ("MMSE_rates", [ :color => "r", :linestyle => "--",  :label => "MaxSINR (MMSE)" ]),
        ],

        "WeightedMaxSINR" => [
            ("logdet_rates", [ :color => "g", :linestyle => "-", :label => "WeightedMaxSINR (logdet)" ]),
            ("MMSE_rates", [ :color => "g", :linestyle => "--",  :label => "WeightedMaxSINR (MMSE)" ]),
        ],

        "Komulainen2013_WMMSE" => [
            ("logdet_rates", [ :color => "y", :linestyle => "-", :label => "DiagWMMSE (logdet)" ]),
            ("MMSE_rates", [ :color => "y", :linestyle => "--",  :label => "DiagWMMSE (MMSE)" ]),
        ],

        "Razaviyayn2013_MinMaxWMMSE" => [
            ("logdet_rates", [ :color => "m", :linestyle => "-", :label => "MinMaxWMMSE (logdet)" ]),
            ("MMSE_rates", [ :color => "m", :linestyle => "--",  :label => "MinMaxWMMSE (MMSE)" ]),
        ],

        "Eigenprecoding" => [
            ("intercell_tdma_logdet_rates", [ :color => "c", :linestyle => "-", :label => "TDMA" ]),
            ("intracell_tdma_logdet_rates", [ :color => "c", :linestyle => "-.",  :label => "Intracell TDMA" ]),
            ("uncoord_logdet_rates", [ :color => "k", :linestyle => "-", :label => "Uncoord. transm." ]),
        ],
    ],
]

##########################################################################
# Plot it
for file_name in parsed_args["file_names"]
    data = load(file_name)
    processed_results = postprocess_convergence(data["raw_results"], data["simulation_params"], plot_params)
    plot_convergence(processed_results, data["simulation_params"], plot_params)
end
