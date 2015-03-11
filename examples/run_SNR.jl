#!/usr/bin/env julia

##########################################################################
# run_SNR.jl
#
# Run all SNR simulations.
##########################################################################
include("run_SNR-ibc.jl")
include("run_SNR-indoors.jl")
include("run_SNR-randomlargescale.jl")
include("run_SNR-triangular3site.jl")
include("run_SNR-triangularhetnet.jl")
