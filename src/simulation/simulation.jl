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
function simulate(network::Network, simulation_params::SimulationParams)
    # Number of drops and small scale fading realizations
    Ndrops = simulation_params["Ndrops"]
    Nsim = simulation_params["Nsim"]

    # Main independent variable
    idp_func = simulation_params["independent_variable"][1]
    idp_vals = simulation_params["independent_variable"][2]
    idp_vals_length = length(idp_vals)

    # Auxiliary independent variables
    Naux = length(simulation_params["aux_independent_variables"])
    aux_idp_funcs = [ simulation_params["aux_independent_variables"][n][1] for n = 1:Naux ]
    aux_idp_vals = [ simulation_params["aux_independent_variables"][n][2] for n = 1:Naux ]

    # Check that all independent variable vectors are the same length
    aux_idp_vals_length = length(simulation_params["aux_independent_variables"][1][2])
    for n = 2:Naux
        aux_idp_vals_length == length(simulation_params["aux_independent_variables"][n][2]) ? nothing : error("Auxiliary independent variable vectors must have equal length.")
    end

    println("-- simulate on $network.")
    println("--- Ndrops: $Ndrops, Nsim: $Nsim.")
    Lumberjack.info("Starting simulation.",
        { :network => network, :simulation_params => simulation_params })

    cell_assignment = assign_cells_by_id(network)

    # Set initial aux precoding params
    set_aux_precoding_params!(network, simulation_params["aux_precoding_params"])

    # Ensure that we are using output protocol 2, so we don't have to store all
    # the intermediate iterations.
    set_aux_precoding_param!(network, 2, "output_protocol")

    # Storage container for results
    raw_results = MultipleSimulationResults(Ndrops, Nsim, idp_vals_length, aux_idp_vals_length)

    # Main simulation loop
    tic()
    for Ndrops_idx = 1:Ndrops
        if Ndrops_idx == Ndrops
            println(Ndrops_idx, ".")
        else
            print(Ndrops_idx, ", ")
        end
        Lumberjack.info("Looping over drop $Ndrops_idx/$Ndrops.")

        draw_user_drop!(network)

        for Nsim_idx = 1:Nsim
            Lumberjack.info("Looping over sim $Nsim_idx/$Nsim.")

            channel = draw_channel(network)

            # Loop over main independent variable
            for idp_vals_idx = 1:idp_vals_length
                # Set main independent variable
                idp_func(network, idp_vals[idp_vals_idx])

                # Loop over auxiliary variables
                for aux_idp_vals_idx = 1:aux_idp_vals_length
                    # Set all auxiliary independent variables
                    for Naux_idx = 1:Naux
                        aux_idp_funcs[Naux_idx](network, aux_idp_vals[Naux_idx][aux_idp_vals_idx])
                    end

                    # Run precoding methods
                    current_results = SingleSimulationResults()
                    for method in simulation_params["precoding_methods"]
                        current_results[string(method)] = method(channel, network, cell_assignment)
                    end
                    raw_results[Ndrops_idx, Nsim_idx, idp_vals_idx, aux_idp_vals_idx] = current_results
                end
            end
        end
    end
    t = toq(); println("--- elapsed time: ", t/60, " minutes")

    return raw_results
end

function process(raw_results::MultipleSimulationResults,
    simulation_params::SimulationParams, plot_params::PlotParams)

    precoding_methods = intersect(simulation_params["precoding_methods"], keys(plot_params["precoding_methods"]))
    Ndrops = simulation_params["Ndrops"]
    Nsim = simulation_params["Nsim"]

    # Length of independent variable vectors
    idp_vals_length = length(simulation_params["independent_variable"][2])
    Naux = length(simulation_params["aux_independent_variables"])
    aux_idp_vals_length = length(simulation_params["aux_independent_variables"][1][2])
    for n = 2:Naux
        aux_idp_vals_length == length(simulation_params["aux_independent_variables"][n][2]) ? nothing : error("Auxiliary independent variable vectors must have equal length.")
    end

    # Compact results into result matrices
    results = [ string(method) => Dict{ASCIIString, Array{Float64}}() for method in simulation_params["precoding_methods"] ]
    for method_name in precoding_methods
        for (result_param,) in plot_params["precoding_methods"][method_name]
            if isa(result_param, ASCIIString)
                result_name = result_param
            else
                (calculator, calculate_from) = result_param
                result_name = string(string(calculator), "_", calculate_from)
            end

            result_dimensions = size(raw_results[1, 1, 1][method_name][result_name])
            result_ranges = [ 1:s for s in result_dimensions ]
            results[method_name][result_name] = Array(Float64, Ndrops, Nsim, idp_vals_length, aux_idp_vals_length, result_dimensions...)
            for Ndrops_idx = 1:Ndrops; for Nsim_idx = 1:Nsim; for idp_vals_idx = 1:idp_vals_length; for aux_idp_vals_idx = 1:aux_idp_vals_length
                if isa(result_param, ASCIIString)
                    result = raw_results[Ndrops_idx, Nsim_idx, idp_vals_idx, aux_idp_vals_idx][method_name][result_name]
                else
                    result = calculator(raw_results[Ndrops_idx, Nsim_idx, idp_vals_idx, aux_idp_vals_idx][method_name][calculate_from])
                end
                results[method_name][result_name][Ndrops_idx, Nsim_idx, idp_vals_idx, aux_idp_vals_idx, result_ranges...] = result
            end; end; end; end
        end
    end

    # Calculate statistics on result matrices
    results_mean = [
        objective_name =>
            [ string(method) => Dict{ASCIIString, Array{Float64}}() for method in simulation_params["precoding_methods"] ]
        for (objective_name,) in plot_params["objectives"]
    ]
    results_var = [
        objective_name =>
            [ string(method) => Dict{ASCIIString, Array{Float64}}() for method in simulation_params["precoding_methods"] ]
        for (objective_name,) in plot_params["objectives"]
    ]
    for (objective_name, (objective_func,)) in plot_params["objectives"]
        for method_name in precoding_methods
            for (result_param,) in plot_params["precoding_methods"][method_name]
                if isa(result_param, ASCIIString)
                    result_name = result_param
                else
                    result_name = string(string(calculator), "_", calculate_from)
                end

                # mean: average over drops and sims
                results_mean[objective_name][method_name][result_name] = squeeze(mean(objective_func(results[method_name][result_name]), 1:2), [1,2,5,6])

                # var: average over sims, estimate var over drops
                results_var[objective_name][method_name][result_name] = squeeze(var(mean(objective_func(results[method_name][result_name]), 2), 1), [1,2,5,6])
            end
        end
    end

    return results, results_mean, results_var
end

function plot(processed_results, simulation_params::SimulationParams,
    plot_params::PlotParams)

    results = processed_results[1]
    results_mean = processed_results[2]
    results_var = processed_results[3]

    ### SYSTEM-LEVEL OBJECTIVE ###
    for (objective_name, (_, objective_params)) in plot_params["objectives"]
        fig, ax = plot_precoding_methods(
                   results_mean[objective_name],
                   results_var[objective_name],
                   simulation_params,
                   plot_params,
                   xvals=simulation_params["independent_variable"][2])

        set_axis_params!(ax, objective_params)

        if displayable("application/pdf")
            display(fig)
        else
            open("SNR_$(simulation_params["name"])_$(plot_params["name_suffix"])_$(objective_name).pdf", "w") do file
                writemime(file, "application/pdf", fig)
            end
        end
        PyPlot.close(fig)
    end
end

##########################################################################
# Performance test
function simulate_performance(network::Network, simulation_params::SimulationParams)
    println("-- performance test on $network.")
    Lumberjack.info("Starting performance test.",
        { :network => network, :simulation_params => simulation_params })

    cell_assignment = assign_cells_by_id(network)

    draw_user_drop!(network)
    channel = draw_channel(network)

    # Make sure things are JITed
    for method in simulation_params["precoding_methods"]
        method(channel, network, cell_assignment)
    end
    println("--- JITing @time macro")
    @time 1

    # Run performance test
    for method in simulation_params["precoding_methods"]
        println("--- Testing performance of ", string(method))
        @time for i = 1:simulation_params["Ntest"]; method(channel, network, cell_assignment); end
    end
end


##########################################################################
# Other functions
include("convergence.jl")
include("visualization.jl")
