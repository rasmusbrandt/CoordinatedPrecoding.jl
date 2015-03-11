#!/usr/bin/env julia

##########################################################################
# run_convergence.jl
#
# Run all convergence simulations.
##########################################################################
include("run_convergence-ibc.jl")
include("run_convergence-indoors.jl")
include("run_convergence-randomlargescale.jl")
include("run_convergence-triangular3site.jl")
include("run_convergence-triangularhetnet.jl")
