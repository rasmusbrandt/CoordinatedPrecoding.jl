function simulate_convergence(network::Network, simulation_params::SimulationParams)
    Ndrops = simulation_params["Ndrops"]
    Nsim = simulation_params["Nsim"]

    println("-- simulate_convergence on $network.")
    println("--- Ndrops: $Ndrops, Nsim: $Nsim.")
    Lumberjack.info("Starting convergence simulation.",
        { :network => network, :simulation_params => simulation_params })

    cell_assignment = assign_cells_by_id(network)

    # Set initial aux precoding params
    set_aux_precoding_params!(network, simulation_params["aux_precoding_params"])

    # Ensure that we are using output protocol 1, because we want to store all
    # intermediate iterations.
    set_aux_precoding_param!(network, 1, "output_protocol")

    raw_results = MultipleSimulationResults(Ndrops, Nsim)

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

            current_results = SingleSimulationResults()
            for method in simulation_params["precoding_methods"]
                current_results[string(method)] = method(channel, network, cell_assignment)
            end

            raw_results[Ndrops_idx, Nsim_idx] = current_results
        end
    end
    t = toq(); println("--- elapsed time: ", t/60, " minutes")

    return raw_results
end

function process_convergence(raw_results::MultipleSimulationResults,
    simulation_params::SimulationParams,
    plot_params::PlotParams)

    precoding_methods = intersect(simulation_params["precoding_methods"], keys(plot_params["precoding_methods"]))
    Ndrops = simulation_params["Ndrops"]
    Nsim = simulation_params["Nsim"]

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
            results[method_name][result_name] = Array(Float64, Ndrops, Nsim, result_dimensions...)
            for Ndrops_idx = 1:Ndrops; for Nsim_idx = 1:Nsim
                if isa(result_param, ASCIIString)
                    result = raw_results[Ndrops_idx, Nsim_idx][method_name][result_name]
                else
                    result = calculator(raw_results[Ndrops_idx, Nsim_idx][method_name][calculate_from])
                end
                results[method_name][result_name][Ndrops_idx, Nsim_idx, result_ranges...] = result
            end; end
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
                results_mean[objective_name][method_name][result_name] = squeeze(mean(objective_func(results[method_name][result_name]), 1:2), 1:4)

                # var: average over sims, estimate var over drops
                results_var[objective_name][method_name][result_name] = squeeze(var(mean(objective_func(results[method_name][result_name]), 2), 1), 1:4)
            end
        end
    end

    return results, results_mean, results_var
end

function plot_convergence(processed_results, simulation_params, plot_params)
    results = processed_results[1]
    results_mean = processed_results[2]
    results_var = processed_results[3]

    ### SYSTEM-LEVEL OBJECTIVE ###
    for (objective_name, (_, objective_params)) in plot_params["objectives"]
        fig, ax = plot_precoding_methods(
                    results_mean[objective_name],
                    results_var[objective_name],
                    simulation_params,
                    plot_params)

        set_axis_params!(ax, objective_params)

        if displayable("application/pdf")
            display(fig)
        else
            open("convergence_$(simulation_params["name"])_$(plot_params["name_suffix"])_$(objective_name).pdf", "w") do file
                writemime(file, "application/pdf", fig)
            end
        end
        PyPlot.close(fig)
    end

    ### USER UTILITIES ###
    K = simulation_params["I"]*simulation_params["Kc"]
    fig = PyPlot.figure(figsize=(6*K,3*length(intersect(simulation_params["precoding_methods"], keys(plot_params["precoding_methods"])))))
    subplot_ind = 1

    for method_name in intersect(simulation_params["precoding_methods"], keys(plot_params["precoding_methods"]))
        for k = 1:K
            ax = fig[:add_subplot](length(intersect(simulation_params["precoding_methods"], keys(plot_params["precoding_methods"]))), K, subplot_ind); subplot_ind += 1

            for (result_param, result_plot_params) in plot_params["precoding_methods"][method_name]
                if isa(result_param, ASCIIString)
                    result = results[method_name][result_param]
                else
                    (calculator, calculate_from) = result_param
                    result = calculator(results[method_name][calculate_from])
                end

                ax[:plot](1:size(result, 5),
                          squeeze(mean(sum(result[:,:,k,:,:], 4), 1:2), 1:4),
                          result_plot_params["key"],
                          label=result_plot_params["legend"])
            end

            if k == 1
                ax[:set_ylabel](method_name)
            end
            if subplot_ind < length(intersect(simulation_params["precoding_methods"], keys(plot_params["precoding_methods"])))
                ax[:set_title]("User $k")
            end
            if subplot_ind-1 > (K-1)*length(intersect(simulation_params["precoding_methods"], keys(plot_params["precoding_methods"])))
                ax[:set_xlabel]("Iteration")
            end
        end
    end

    if displayable("application/pdf")
        display(fig)
    else
        open("convergence_$(simulation_params["name"])_$(plot_params["name_suffix"])_peruser.pdf", "w") do file
            writemime(file, "application/pdf", fig)
        end
    end
    PyPlot.close(fig)

    ### STREAM UTILITIES ###
    fig = PyPlot.figure(figsize=(6*simulation_params["d"],3*K))
    subplot_ind = 1

    for k = 1:K
        for n = 1:simulation_params["d"]
            ax = fig[:add_subplot](K, simulation_params["d"], subplot_ind); subplot_ind += 1

            for method_name in intersect(simulation_params["precoding_methods"], keys(plot_params["precoding_methods"]))
                for (result_param, result_plot_params) in plot_params["precoding_methods"][method_name]
                    if isa(result_param, ASCIIString)
                        result = results[method_name][result_param]
                    else
                        (calculator, calculate_from) = result_param
                        result = calculator(results[method_name][calculate_from])
                    end

                    ax[:plot](1:size(result, 5),
                              squeeze(mean(result[:,:,k,n,:], 1:2), 1:4),
                              result_plot_params["key"],
                              label=result_plot_params["legend"])
                end
            end

            if n == 1
                ax[:set_ylabel]("User $k")
            end
            if k == 1
                ax[:set_title]("Stream $n")
            end
            if k == K
                ax[:set_xlabel]("Iteration")
            end
        end
    end

    if displayable("application/pdf")
        display(fig)
    else
        open("convergence_$(simulation_params["name"])_$(plot_params["name_suffix"])_perstream.pdf", "w") do file
            writemime(file, "application/pdf", fig)
        end
    end
    PyPlot.close(fig)
end
