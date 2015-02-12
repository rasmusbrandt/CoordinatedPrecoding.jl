function plot_precoding_methods(results_mean, results_var, simulation_params,
    plot_params; xvals=[])

    fig = PyPlot.figure(figsize=plot_params["figsize"])
    ax = fig[:add_subplot](1, 1, 1)

    # Set type of scales for the axes
    haskey(plot_params, "xscale") ? ax[:set_xscale](plot_params["xscale"]) : nothing
    haskey(plot_params, "yscale") ? ax[:set_yscale](plot_params["yscale"]) : nothing

    for method_name in intersect(simulation_params["precoding_methods"], keys(plot_params["precoding_methods"]))
        for (result_param, result_plot_params) in plot_params["precoding_methods"][method_name]
            if isa(result_param, ASCIIString)
                result_name = result_param
            else
                (calculator, calculate_from) = result_param
                result_name = string(string(calculator), "_", calculate_from)
            end

            xvals_plot = (xvals != []) ? xvals : 1:size(results_mean[method_name][result_name], 1)

            ax[:plot](xvals_plot,
                      results_mean[method_name][result_name],
                      result_plot_params["key"],
                      label=result_plot_params["legend"])

            if haskey(plot_params, "confidence_interval_z_alpha_half")
                facecolor = match(r"[a-z]", result_plot_params["key"]).match
                ax[:fill_between](xvals_plot,
                                  results_mean[method_name][result_name] + plot_params["confidence_interval_z_alpha_half"]*results_var[method_name][result_name]/sqrt(simulation_params["Ndrops"]),
                                  results_mean[method_name][result_name] - plot_params["confidence_interval_z_alpha_half"]*results_var[method_name][result_name]/sqrt(simulation_params["Ndrops"]),
                                  facecolor=facecolor, alpha=0.5)
            end
        end
    end

    return (fig, ax)
end

function set_axis_params!(axis, objective_params)
    haskey(objective_params, "xlabel") ? axis[:set_xlabel](objective_params["xlabel"]) : nothing
    haskey(objective_params, "ylabel") ? axis[:set_ylabel](objective_params["ylabel"]) : nothing

    if haskey(objective_params, "legend_loc")
        if haskey(objective_params, "legend_prop")
            axis[:legend](loc=objective_params["legend_loc"],
                        prop=objective_params["legend_prop"])
        else
            axis[:legend](loc=objective_params["legend_loc"])
        end
    end
end

##########################################################################
# General simulation plot
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
            open("$(simulation_params["name"])_$(plot_params["name_suffix"])_$(objective_name).pdf", "w") do file
                writemime(file, "application/pdf", fig)
            end
        end
        PyPlot.close(fig)
    end
end

##########################################################################
# Convergence simulation plot
function plot_convergence(processed_results,
    simulation_params::SimulationParams, plot_params::PlotParams)

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
            open("$(simulation_params["name"])_$(plot_params["name_suffix"])_$(objective_name).pdf", "w") do file
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

                ax[:plot](1:size(result, 6),
                          transpose(squeeze(mean(sum(result[:,:,:,k,:,:], 5), 1:2), [1,2,4,5])),
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
        open("$(simulation_params["name"])_$(plot_params["name_suffix"])_peruser.pdf", "w") do file
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

                    ax[:plot](1:size(result, 6),
                              transpose(squeeze(mean(result[:,:,:,k,n,:], 1:2), [1,2,4,5])),
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
        open("$(simulation_params["name"])_$(plot_params["name_suffix"])_perstream.pdf", "w") do file
            writemime(file, "application/pdf", fig)
        end
    end
    PyPlot.close(fig)
end
