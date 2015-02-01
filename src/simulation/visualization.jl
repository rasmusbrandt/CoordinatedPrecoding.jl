function plot_precoding_methods(results_mean, results_var, simulation_params,
    plot_params; xvals=[])

    fig = PyPlot.figure(figsize=plot_params["figsize"])
    ax = fig[:add_subplot](1, 1, 1)

    for method_name in intersect(simulation_params["precoding_methods"], keys(plot_params["precoding_methods"]))
        for (result_param, result_plot_params) in plot_params["precoding_methods"][method_name]
            if isa(result_param, ASCIIString)
                result_name = result_param
            else
                (calculator, calculate_from) = result_param
                result_name = string(string(calculator), "_", calculate_from)
            end

            xvals_plot = (xvals != []) ? xvals : 1:length(results_mean[method_name][result_name])

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
