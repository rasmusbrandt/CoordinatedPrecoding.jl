##########################################################################
# Type for simulation parameters
typealias SimulationParams Dict{ASCIIString, Any}

##########################################################################
# Simulation functions
include("convergence.jl")
include("performance.jl")
include("SNR.jl")
