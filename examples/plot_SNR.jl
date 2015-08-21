#!/usr/bin/env julia

##########################################################################
# plot_SNR.jl
#
# Plots SNR curves.
##########################################################################

using CoordinatedPrecoding
using Compat, JLD, ArgParse

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
plot_params = @compat Dict(
    "plot_name" => "",

    "objective" => :sumrate,

    "figure" => @Compat.Dict(
        :figsize => (8,5),
        :dpi => 125,
        :facecolor => "w",
        :edgecolor => "k",
    ),

    "axes" => @Compat.Dict(
        :xlabel => "Transmit power [dBm]",
        :ylabel => "Sum rate [bits/s/Hz]",
        :xscale => "linear",
        :yscale => "linear",
    ),

    "legend" => @Compat.Dict(
        :loc => "best",
        :fontsize => 8,
    ),

    # "confidence_interval_zalpha2" => 1.96,

    "methods" => @Compat.Dict(
        "Shi2011_WMMSE" => [
            ("logdet_rates", @compat Dict(:color => "b", :linestyle => "-", :label => "WMMSE (logdet)")),
            ("MMSE_rates", @compat Dict(:color => "b", :linestyle => "--",  :label => "WMMSE (MMSE)")),
        ],

        "Gomadam2008_MaxSINR" => [
            ("logdet_rates", @compat Dict(:color => "r", :linestyle => "-", :label => "MaxSINR (logdet)")),
            ("MMSE_rates", @compat Dict(:color => "r", :linestyle => "--",  :label => "MaxSINR (MMSE)")),
        ],

        "Gomadam2008_MinWLI" => [
            ("logdet_rates", @compat Dict(:color => "gray", :linestyle => "-", :label => "MinWLI (logdet)")),
            ("MMSE_rates", @compat Dict(:color => "gray", :linestyle => "--",  :label => "MinWLI (MMSE)")),
        ],

        "Komulainen2013_WMMSE" => [
            ("logdet_rates", @compat Dict(:color => "y", :linestyle => "-", :label => "DiagWMMSE (logdet)")),
            ("MMSE_rates", @compat Dict(:color => "y", :linestyle => "--",  :label => "DiagWMMSE (MMSE)")),
        ],

        "Razaviyayn2013_MinMaxWMMSE" => [
            ("logdet_rates", @compat Dict(:color => "m", :linestyle => "-", :label => "MinMaxWMMSE (logdet)")),
            ("MMSE_rates", @compat Dict(:color => "m", :linestyle => "--",  :label => "MinMaxWMMSE (MMSE)")),
        ],

        "Eigenprecoding" => [
            ("intercell_tdma_logdet_rates", @compat Dict(:color => "c", :linestyle => "-", :label => "TDMA")),
            ("intracell_tdma_logdet_rates", @compat Dict(:color => "c", :linestyle => "-.",  :label => "Intracell TDMA")),
            ("uncoord_logdet_rates", @compat Dict(:color => "k", :linestyle => "-", :label => "Uncoord. transm.")),
        ],
    ),
)

##########################################################################
# Plot it
for file_name in parsed_args["file_names"]
    data = load(file_name)
    processed_results = postprocess(data["raw_results"], data["simulation_params"], plot_params)
    plot(processed_results, data["simulation_params"], plot_params)
end
