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
# General simulation functions
function simulate(network::Network, simulation_params::SimulationParams; loop_over = :precoding_methods)
    # Number of drops and small scale fading realizations
    Ndrops = simulation_params["Ndrops"]
    Nsim = simulation_params["Nsim"]

    # Main independent variable
    idp_func = simulation_params["independent_variable"][1]
    idp_vals = simulation_params["independent_variable"][2]
    idp_vals_length = length(idp_vals)

    # Auxiliary independent variables
    Naux, aux_idp_funcs, aux_idp_vals, aux_idp_vals_length = get_aux_idp(simulation_params)

    # Get assignment/precoding functions depending on looping mode
    precoding_method, assignment_method = get_other_method(simulation_params, loop_over)

    println("-- simulate on $network.")
    println("--- Ndrops: $Ndrops, Nsim: $Nsim.")
    Lumberjack.info("Starting simulation.",
        { :network => network, :simulation_params => simulation_params })

    # Set initial aux params
    set_initial_aux_params!(simulation_params, network)

    # Ensure that we are not storing all intermediate iterations.
    set_aux_precoding_param!(network, :final_iteration, "output_protocol")

    # Storage container for results
    raw_results = MultipleSimulationResults(Ndrops, Nsim, idp_vals_length, aux_idp_vals_length)

    # Outer simulation loop
    if loop_over == :precoding_methods
        progress = ProgressMeter.Progress(Ndrops*idp_vals_length*aux_idp_vals_length*length(simulation_params["precoding_methods"])*Nsim)
    else
        progress = ProgressMeter.Progress(Ndrops*idp_vals_length*aux_idp_vals_length*length(simulation_params["assignment_methods"])*Nsim)
    end
    for Ndrops_idx = 1:Ndrops
        draw_user_drop!(network)
        channels = [ draw_channel(network) for n = 1:Nsim ]

        # Inner simulation loop depends on loop mode
        if loop_over == :precoding_methods
            # Loop over main independent variable
            for idp_vals_idx = 1:idp_vals_length
                # Set main independent variable
                idp_func(network, idp_vals[idp_vals_idx])

                # Loop over auxiliary variables
                for aux_idp_vals_idx = 1:aux_idp_vals_length
                    if Naux != 0
                        # Set all auxiliary independent variables
                        for Naux_idx = 1:Naux
                            aux_idp_funcs[Naux_idx](network, aux_idp_vals[Naux_idx][aux_idp_vals_idx])
                        end
                    end

                    assignment_method(channels[1], network)

                    for precoding_method in simulation_params["precoding_methods"]
                        for Nsim_idx = 1:Nsim
                            # Allocate memory if this is the first method to be run
                            if precoding_method == simulation_params["precoding_methods"][1]
                                raw_results[Ndrops_idx, Nsim_idx, idp_vals_idx, aux_idp_vals_idx] = SingleSimulationResults()
                            end
                            raw_results[Ndrops_idx, Nsim_idx, idp_vals_idx, aux_idp_vals_idx][string(precoding_method)] = precoding_method(channels[Nsim_idx], network)

                            ProgressMeter.next!(progress)
                        end
                    end
                end
            end
        elseif loop_over == :assignment_methods
            # Loop over main independent variable
            for idp_vals_idx = 1:idp_vals_length
                # Set main independent variable
                idp_func(network, idp_vals[idp_vals_idx])

                # Loop over auxiliary variables
                for aux_idp_vals_idx = 1:aux_idp_vals_length
                    if Naux != 0
                        # Set all auxiliary independent variables
                        for Naux_idx = 1:Naux
                            aux_idp_funcs[Naux_idx](network, aux_idp_vals[Naux_idx][aux_idp_vals_idx])
                        end
                    end

                    for assignment_method in simulation_params["assignment_methods"]
                        assignment_method(channels[1], network)

                        for Nsim_idx = 1:Nsim
                            # Allocate memory if this is the first method to be run
                            if assignment_method == simulation_params["assignment_methods"][1]
                                raw_results[Ndrops_idx, Nsim_idx, idp_vals_idx, aux_idp_vals_idx] = SingleSimulationResults()
                            end
                            raw_results[Ndrops_idx, Nsim_idx, idp_vals_idx, aux_idp_vals_idx][string(assignment_method)] = precoding_method(channels[Nsim_idx], network)

                            ProgressMeter.next!(progress)
                        end
                    end
                end
            end
        end
    end

    return raw_results
end

##########################################################################
# Convergence simulation functions
function simulate_convergence(network::Network, simulation_params::SimulationParams; loop_over = :precoding_methods)
    # Number of drops and small scale fading realizations
    Ndrops = simulation_params["Ndrops"]
    Nsim = simulation_params["Nsim"]

    # Auxiliary independent variables
    Naux, aux_idp_funcs, aux_idp_vals, aux_idp_vals_length = get_aux_idp(simulation_params)

    # Get assignment/precoding functions depending on looping mode
    precoding_method, assignment_method = get_other_method(simulation_params, loop_over)

    println("-- simulate_convergence on $network.")
    println("--- Ndrops: $Ndrops, Nsim: $Nsim.")
    Lumberjack.info("Starting convergence simulation.",
        { :network => network, :simulation_params => simulation_params })

    # Set initial aux params
    set_initial_aux_params!(simulation_params, network)

    # We want to store all intermediate iterations.
    set_aux_precoding_param!(network, :all_iterations, "output_protocol")

    # Warn if we are not using a fixed number of iterations.
    if get_aux_precoding_param(network, "stop_crit") > 0.
        Lumberjack.warn("process_convergence will not be able to run since stop_crit is non-zero, and the algorithms may therefore use different numbers of iterations.")
    end

    # Storage container for results
    raw_results = MultipleSimulationResults(Ndrops, Nsim, aux_idp_vals_length)

    # Outer simulation loop
    progress = ProgressMeter.Progress(Ndrops*aux_idp_vals_length*length(simulation_params["precoding_methods"])*Nsim)
    for Ndrops_idx = 1:Ndrops
        draw_user_drop!(network)
        channels = [ draw_channel(network) for n = 1:Nsim ]

        # Inner simulation loop depends on loop mode
        if loop_over == :precoding_methods
            # Loop over auxiliary variables
            for aux_idp_vals_idx = 1:aux_idp_vals_length
                if Naux != 0
                    # Set all auxiliary independent variables
                    for Naux_idx = 1:Naux
                        aux_idp_funcs[Naux_idx](network, aux_idp_vals[Naux_idx][aux_idp_vals_idx])
                    end
                end

                assignment_method(channels[1], network)

                for precoding_method in simulation_params["precoding_methods"]
                    for Nsim_idx = 1:Nsim
                        # Allocate memory if this is the first method to be run
                        if precoding_method == simulation_params["precoding_methods"][1]
                            raw_results[Ndrops_idx, Nsim_idx, aux_idp_vals_idx] = SingleSimulationResults()
                        end
                        raw_results[Ndrops_idx, Nsim_idx, aux_idp_vals_idx][string(precoding_method)] = precoding_method(channels[Nsim_idx], network)

                        ProgressMeter.next!(progress)
                    end
                end
            end
        elseif loop_over == :assignment_methods
            # Loop over auxiliary variables
            for aux_idp_vals_idx = 1:aux_idp_vals_length
                if Naux != 0
                    # Set all auxiliary independent variables
                    for Naux_idx = 1:Naux
                        aux_idp_funcs[Naux_idx](network, aux_idp_vals[Naux_idx][aux_idp_vals_idx])
                    end
                end

                for assignment_method in simulation_params["assignment_methods"]
                    assignment_method(channels[1], network)

                    for Nsim_idx = 1:Nsim
                        # Allocate memory if this is the first method to be run
                        if assignment_method == simulation_params["assignment_methods"][1]
                            raw_results[Ndrops_idx, Nsim_idx, aux_idp_vals_idx] = SingleSimulationResults()
                        end
                        raw_results[Ndrops_idx, Nsim_idx, aux_idp_vals_idx][string(assignment_method)] = precoding_method(channels[Nsim_idx], network)

                        ProgressMeter.next!(progress)
                    end
                end
            end
        end
    end

    return raw_results
end

##########################################################################
# Performance test
function simulate_performance(network::Network, simulation_params::SimulationParams)
    println("-- performance test on $network.")
    Lumberjack.info("Starting performance test.",
        { :network => network, :simulation_params => simulation_params })

    assign_cells_by_id!(network)

    # Set initial aux precoding params
    set_aux_precoding_params!(network, simulation_params["aux_precoding_params"])

    # No point storing all intermediate iterations.
    set_aux_precoding_param!(network, :final_iteration, "output_protocol")

    draw_user_drop!(network)
    channel = draw_channel(network)

    # Make sure things are JITed
    for method in simulation_params["precoding_methods"]
        method(channel, network)
    end
    println("--- JITing @time macro")
    @time 1

    # Run performance test
    for method in simulation_params["precoding_methods"]
        println("--- Testing performance of ", string(method))
        @time for i = 1:simulation_params["Ntest"]; method(channel, network); end
    end
end

function get_other_method(simulation_params, loop_over)
    # Dummies
    precoding_method() = nothing
    assignment_method() = nothing

    if loop_over == :precoding_methods
        if haskey(simulation_params, "assignment_methods")
            if length(simulation_params["assignment_methods"]) > 1
                warn("Looping over precoding methods: will only use first assignment method provided.")
            end
            assignment_method(channel, network) = simulation_params["assignment_methods"][1](channel, network)
        else
            assignment_method(channel, network) = assign_cells_by_id!(network)
        end
    elseif loop_over == :assignment_methods
        if haskey(simulation_params, "precoding_methods")
            if length(simulation_params["precoding_methods"]) > 1
                warn("Looping over assignment methods: will only use first precoding method provided.")
            end
            precoding_method(channel, network) = simulation_params["precoding_methods"][1](channel, network)
        else
            error("No precoding method specified.")
        end
    end

    return precoding_method, assignment_method
end

function get_aux_idp(simulation_params)
    if haskey(simulation_params, "aux_independent_variables")
        Naux = length(simulation_params["aux_independent_variables"])
        aux_idp_funcs = [ simulation_params["aux_independent_variables"][n][1] for n = 1:Naux ]
        aux_idp_vals = [ simulation_params["aux_independent_variables"][n][2] for n = 1:Naux ]

        # Check that all independent variable vectors are the same length
        aux_idp_vals_length = length(simulation_params["aux_independent_variables"][1][2])
        for n = 2:Naux
            aux_idp_vals_length == length(simulation_params["aux_independent_variables"][n][2]) ? nothing : error("Auxiliary independent variable vectors must have equal length.")
        end
    else
        Naux = 0
        aux_idp_funcs = []
        aux_idp_vals = []
        aux_idp_vals_length = 1
    end

    return Naux, aux_idp_funcs, aux_idp_vals, aux_idp_vals_length
end

function set_initial_aux_params!(simulation_params, network)
    haskey(simulation_params, "aux_precoding_params") && set_aux_precoding_params!(network, simulation_params["aux_precoding_params"])
    haskey(simulation_params, "aux_assignment_params") && set_aux_assignment_params!(network, simulation_params["aux_assignment_params"])
end

##########################################################################
# Other functions
include("visualization.jl")
