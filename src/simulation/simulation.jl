##########################################################################
# Simulation types
typealias SimulationParams Dict{ASCIIString, Any}
typealias PlotParams Dict{ASCIIString, Any}

type SingleSimulationResults
    precoding_results::Dict{ASCIIString, PrecodingResults}
end
SingleSimulationResults() =
    SingleSimulationResults(Dict{ASCIIString, PrecodingResults}())
Base.getindex(s::SingleSimulationResults, k::ASCIIString) =
    getindex(s.precoding_results, k)
Base.setindex!(s::SingleSimulationResults, v, k::ASCIIString) =
    setindex!(s.precoding_results, v, k)

type MultipleSimulationResults{N}
    simulation_results::Array{SingleSimulationResults, N}
end
MultipleSimulationResults(dims...) =
    MultipleSimulationResults(Array(SingleSimulationResults, dims...))
Base.getindex(m::MultipleSimulationResults, inds...) =
    getindex(m.simulation_results, inds...)
Base.setindex!(m::MultipleSimulationResults, v::SingleSimulationResults, inds...) =
    setindex!(m.simulation_results, v, inds...)

##########################################################################
# Simulation functions
include("convergence.jl")
include("performance.jl")
include("SNR.jl")
include("visualization.jl")
