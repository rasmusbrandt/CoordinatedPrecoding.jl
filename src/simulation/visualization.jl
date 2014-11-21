function plot_methods(results_mean, results_var, simulation_params,
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

function set_axis_params!(axis, systemlevel_params)
    haskey(systemlevel_params, "xlabel") ? axis[:set_xlabel](systemlevel_params["xlabel"]) : nothing
    haskey(systemlevel_params, "ylabel") ? axis[:set_ylabel](systemlevel_params["ylabel"]) : nothing

    if haskey(systemlevel_params, "legend_loc")
        if haskey(systemlevel_params, "legend_prop")
            axis[:legend](loc=systemlevel_params["legend_loc"],
                        prop=systemlevel_params["legend_prop"])
        else
            axis[:legend](loc=systemlevel_params["legend_loc"])
        end
    end
end

function plot_convergence(results, simulation_params, precoding_settings, plot_params)
    ### SYSTEM-LEVEL POSTPROCESSING ###
    results_mean = [
        systemlevel_name =>
            [ string(method) => Dict{ASCIIString, Array{Float64}}() for method in simulation_params["precoding_methods"] ]
        for (systemlevel_name,) in plot_params["systemlevel_objectives"]
    ]
    results_var = [
        systemlevel_name =>
            [ string(method) => Dict{ASCIIString, Array{Float64}}() for method in simulation_params["precoding_methods"] ]
        for (systemlevel_name,) in plot_params["systemlevel_objectives"]
    ]

    for (systemlevel_name, (systemlevel_func,)) in plot_params["systemlevel_objectives"]
        for method_name in intersect(simulation_params["precoding_methods"], keys(plot_params["precoding_methods"]))
            for (result_param,) in plot_params["precoding_methods"][method_name]
                if isa(result_param, ASCIIString)
                    result_name = result_param
                    result = results[method_name][result_name]
                else
                    (calculator, calculate_from) = result_param
                    result_name = string(string(calculator), "_", calculate_from)
                    result = calculator(results[method_name][calculate_from])
                end

                # mean: average over drops and sims
                # var: average over sims, estimate var over drops
                results_mean[systemlevel_name][method_name][result_name] = squeeze(mean(systemlevel_func(result), 1:2), 1:4)
                results_var[systemlevel_name][method_name][result_name] = squeeze(var(mean(systemlevel_func(result), 2), 1), 1:4)
            end
        end
    end

    ### SYSTEM-LEVEL OBJECTIVE EVOLUTION ###
    for (systemlevel_name, (_, systemlevel_params)) in plot_params["systemlevel_objectives"]
        fig, ax = plot_methods(results_mean[systemlevel_name],
                               results_var[systemlevel_name],
                               simulation_params,
                               plot_params)

        set_axis_params!(ax, systemlevel_params)

        if displayable("application/pdf")
            display(fig)
        else
            open("convergence_$(simulation_params["name"])_$(systemlevel_name).pdf", "w") do file
                writemime(file, "application/pdf", fig)
            end
        end
        PyPlot.close(fig)
    end

    ### USER EVOLUTION ###
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
            if subplot_ind-1 < length(intersect(simulation_params["precoding_methods"], keys(plot_params["precoding_methods"])))
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
        open("convergence_$(simulation_params["name"])_peruser.pdf", "w") do file
            writemime(file, "application/pdf", fig)
        end
    end
    PyPlot.close(fig)

    ### STREAM EVOLUTION ###
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
        open("convergence_$(simulation_params["name"])_perstream.pdf", "w") do file
            writemime(file, "application/pdf", fig)
        end
    end
    PyPlot.close(fig)
end

function plot_SNR(results, simulation_params, precoding_settings, plot_params)
    ### SYSTEM-LEVEL POSTPROCESSING ###
    results_mean = [
        systemlevel_name =>
            [ string(method) => Dict{ASCIIString, Array{Float64}}() for method in simulation_params["precoding_methods"] ]
        for (systemlevel_name,) in plot_params["systemlevel_objectives"]
    ]
    results_var = [
        systemlevel_name =>
            [ string(method) => Dict{ASCIIString, Array{Float64}}() for method in simulation_params["precoding_methods"] ]
        for (systemlevel_name,) in plot_params["systemlevel_objectives"]
    ]

    for (systemlevel_name, (systemlevel_func,)) in plot_params["systemlevel_objectives"]
        for method_name in intersect(simulation_params["precoding_methods"], keys(plot_params["precoding_methods"]))
            for (result_param,) in plot_params["precoding_methods"][method_name]
                if isa(result_param, ASCIIString)
                    result_name = result_param
                    result = results[method_name][result_name]
                else
                    (calculator, calculate_from) = result_param
                    result_name = string(string(calculator), "_", calculate_from)
                    result = calculator(results[method_name][calculate_from])
                end

                # mean: average over drops and sims
                # var: average over sims, estimate var over drops
                results_mean[systemlevel_name][method_name][result_name] = squeeze(mean(systemlevel_func(result), 1:2), [1,2,4,5])
                results_var[systemlevel_name][method_name][result_name] = squeeze(var(mean(systemlevel_func(result), 2), 1), [1,2,4,5])
            end
        end
    end

    ### SYSTEM-LEVEL OBJECTIVE EVOLUTION ###
    for (systemlevel_name, (_, systemlevel_params)) in plot_params["systemlevel_objectives"]
        fig, ax = plot_methods(results_mean[systemlevel_name],
                               results_var[systemlevel_name],
                               simulation_params,
                               plot_params,
                               xvals=simulation_params["Ps_dBm"])

        set_axis_params!(ax, systemlevel_params)

        if displayable("application/pdf")
            display(fig)
        else
            open("SNR_$(simulation_params["name"])_$(systemlevel_name).pdf", "w") do file
                writemime(file, "application/pdf", fig)
            end
        end
        PyPlot.close(fig)
    end
end
