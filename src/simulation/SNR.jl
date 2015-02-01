function simulate_SNR(network::Network,
    simulation_params::SimulationParams,
    precoding_params::PrecodingParams)

    Ndrops = simulation_params["Ndrops"]
    Nsim = simulation_params["Nsim"]
    Npowers = length(simulation_params["Ps_dBm"])

    println("-- simulate_SNR on $network.")
    println("--- Ndrops: $Ndrops, Nsim: $Nsim.")
    Lumberjack.info("Starting SNR simulation.",
        { :network => network, :simulation_params => simulation_params,
          :precoding_params => precoding_params })

    cell_assignment = assign_cells_by_id(network)

    # Ensure that we are using output protocol 2, so we don't have to store all
    # the intermediate iterations.
    precoding_params["output_protocol"] = 2

    raw_results = MultipleSimulationResults(Ndrops, Nsim, Npowers)

    tic()
    for Ndrops_idx = 1:Ndrops
        if Ndrops_idx == Ndrops
            println(Ndrops_idx, ".")
        else
            print(Ndrops_idx, ", ")
        end

        draw_user_drop!(network)

        for Nsim_idx = 1:Nsim
            channel = draw_channel(network)

            for Npowers_idx = 1:Npowers
                set_transmit_powers!(network, 10^(simulation_params["Ps_dBm"][Npowers_idx]/10))

                current_results = SingleSimulationResults()
                for method in simulation_params["precoding_methods"]
                    current_results[string(method)] = method(channel, network, cell_assignment, precoding_params)
                end

                raw_results[Ndrops_idx, Nsim_idx, Npowers_idx] = current_results
            end
        end
    end
    t = toq(); println("--- elapsed time: ", t/60, " minutes")

    return raw_results
end

function process_SNR(raw_results::MultipleSimulationResults,
    simulation_params::SimulationParams,
    plot_params::PlotParams)

    precoding_methods = intersect(simulation_params["precoding_methods"], keys(plot_params["precoding_methods"]))
    Ndrops = simulation_params["Ndrops"]
    Nsim = simulation_params["Nsim"]
    Npowers = length(simulation_params["Ps_dBm"])

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
            results[method_name][result_name] = Array(Float64, Ndrops, Nsim, Npowers, result_dimensions...)
            for Ndrops_idx = 1:Ndrops; for Nsim_idx = 1:Nsim; for Npowers_idx = 1:Npowers
                if isa(result_param, ASCIIString)
                    result = raw_results[Ndrops_idx, Nsim_idx, Npowers_idx][method_name][result_name]
                else
                    result = calculator(raw_results[Ndrops_idx, Nsim_idx, Npowers_idx][method_name][calculate_from])
                end
                results[method_name][result_name][Ndrops_idx, Nsim_idx, Npowers_idx, result_ranges...] = result
            end; end; end
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
                results_mean[objective_name][method_name][result_name] = squeeze(mean(objective_func(results[method_name][result_name]), 1:2), [1,2,4,5])

                # var: average over sims, estimate var over drops
                results_var[objective_name][method_name][result_name] = squeeze(var(mean(objective_func(results[method_name][result_name]), 2), 1), [1,2,4,5])
            end
        end
    end

    return results, results_mean, results_var
end

function plot_SNR(processed_results, simulation_params, plot_params)
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
                   xvals=simulation_params["Ps_dBm"])

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
