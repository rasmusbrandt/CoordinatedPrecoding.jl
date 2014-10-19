function plot_convergence(results, simulation_params, precoding_settings,
    plot_params; plot_confidence_interval = false, z_alpha_half = 1.96)

    ### POSTPROCESSING ###
    results_sumrate_mean = [ string(method) => Dict{ASCIIString, Array{Float64}}() for method in simulation_params["precoding_methods"] ]
    results_sumrate_var  = [ string(method) => Dict{ASCIIString, Array{Float64}}() for method in simulation_params["precoding_methods"] ]
    results_minrate_mean = [ string(method) => Dict{ASCIIString, Array{Float64}}() for method in simulation_params["precoding_methods"] ]
    results_minrate_var  = [ string(method) => Dict{ASCIIString, Array{Float64}}() for method in simulation_params["precoding_methods"] ]

    for method_name in intersect(simulation_params["precoding_methods"], keys(plot_params["precoding_methods"]))
        for (result_param, plot_key, plot_legend) in plot_params["precoding_methods"][method_name]
            if isa(result_param, ASCIIString)
                result_name = result_param
                result = results[method_name][result_name]
            else
                (calculator, calculate_from) = result_param
                result_name = string(string(calculator), "_", calculate_from)
                result = calculator(results[method_name][calculate_from])
            end

            if ndims(result) == 5
                # Ndrops-by-Nsim-by-K-by-d-by-stop_crit

                # mean: average over drops and sims
                results_sumrate_mean[method_name][result_name] = squeeze(mean(sum(result, 3:4), 1:2), 1:4)
                results_minrate_mean[method_name][result_name] = squeeze(mean(minimum(sum(result, 4), 3), 1:2), 1:4)

                # var: average over sims, estimate var over drops
                results_sumrate_var[method_name][result_name] = squeeze(var(mean(sum(result, 3:4), 2), 1), 1:4)
                results_minrate_var[method_name][result_name] = squeeze(var(mean(minimum(sum(result, 4), 3), 2), 1), 1:4)
            elseif ndims(result) == 3
                # Ndrops-by-Nsim-by-stop_crit

                # mean: average over drops and sims
                results_sumrate_mean[method_name][result_name] = squeeze(mean(result, 1:2), 1:2)
                results_minrate_mean[method_name][result_name] = squeeze(mean(result, 1:2), 1:2)

                # var: average over sims, estimate var over drops
                results_sumrate_var[method_name][result_name] = squeeze(var(mean(result, 2), 1), 1:2)
                results_minrate_var[method_name][result_name] = squeeze(var(mean(result, 2), 1), 1:2)
            end
        end
    end

    ### SUM RATE EVOLUTION ###
    fig, ax = plot_methods(1:precoding_settings["stop_crit"],
                results_sumrate_mean, results_sumrate_var, simulation_params,
                plot_params, plot_confidence_interval, z_alpha_half)

    #ax[:legend](loc="lower right")
    ax[:set_xlabel]("Iteration")
    ax[:set_ylabel]("Sum rate [bits/s/Hz]")

    if displayable("application/pdf")
        display(fig)
    else
        open("convergence_$(simulation_params["name"])_sumrate.pdf", "w") do file
            writemime(file, "application/pdf", fig)
        end
    end
    PyPlot.close(fig)

    ### MIN RATE EVOLUTION ###
    fig, ax = plot_methods(1:precoding_settings["stop_crit"],
                results_minrate_mean, results_minrate_var, simulation_params,
                plot_params, plot_confidence_interval, z_alpha_half)

    #ax[:legend](loc="lower right")
    ax[:set_xlabel]("Iteration")
    ax[:set_ylabel]("Min rate [bits/s/Hz]")

    if displayable("application/pdf")
        display(fig)
    else
        open("convergence_$(simulation_params["name"])_minrate.pdf", "w") do file
            writemime(file, "application/pdf", fig)
        end
    end
    PyPlot.close(fig)

    ### USER RATE EVOLUTION ###
    K = simulation_params["I"]*simulation_params["Kc"]
    fig = PyPlot.figure(figsize=(6*K,3*length(intersect(simulation_params["precoding_methods"], keys(plot_params["precoding_methods"])))))
    subplot_ind = 1

    title_printed = false
    for method_name in intersect(simulation_params["precoding_methods"], keys(plot_params["precoding_methods"]))
        for k = 1:K
            ax = fig[:add_subplot](length(intersect(simulation_params["precoding_methods"], keys(plot_params["precoding_methods"]))), K, subplot_ind); subplot_ind += 1

            for (result_param, plot_key, plot_legend) in plot_params["precoding_methods"][method_name]
                if isa(result_param, ASCIIString)
                    result = results[method_name][result_param]
                else
                    (calculator, calculate_from) = result_param
                    result = calculator(results[method_name][calculate_from])
                end

                if ndims(result) == 5
                    # Ndrops-by-Nsim-by-K-by-d-by-stop_crit
                    ax[:plot](1:precoding_settings["stop_crit"], squeeze(mean(sum(result[:,:,k,:,:], 4), 1:2), 1:4), plot_key, label=plot_legend)
                end
            end

            if k == 1
                ax[:set_ylabel](method_name)
            end
            if !title_printed
                ax[:set_title]("User $k")
            end
            if subplot_ind-1 >= (K-1)*length(intersect(simulation_params["precoding_methods"], keys(plot_params["precoding_methods"])))
                ax[:set_xlabel]("Iteration")
            end
        end
        title_printed = true
    end

    if displayable("application/pdf")
        display(fig)
    else
        open("convergence_$(simulation_params["name"])_userrate.pdf", "w") do file
            writemime(file, "application/pdf", fig)
        end
    end
    PyPlot.close(fig)

    ### STREAM RATE EVOLUTION ###
    fig = PyPlot.figure(figsize=(6*simulation_params["d"],3*K))
    subplot_ind = 1

    for k = 1:K
        for n = 1:simulation_params["d"]
            ax = fig[:add_subplot](K, simulation_params["d"], subplot_ind); subplot_ind += 1

            for method_name in intersect(simulation_params["precoding_methods"], keys(plot_params["precoding_methods"]))
                for (result_param, plot_key, plot_legend) in plot_params["precoding_methods"][method_name]
                    if isa(result_param, ASCIIString)
                        result = results[method_name][result_param]
                    else
                        (calculator, calculate_from) = result_param
                        result = calculator(results[method_name][calculate_from])
                    end

                    if ndims(result) == 5
                        # Ndrops-by-Nsim-by-K-by-d-by-stop_crit
                        ax[:plot](1:precoding_settings["stop_crit"], squeeze(mean(result[:,:,k,n,:], 1:2), 1:4), plot_key, label=plot_legend)
                    end
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
        open("convergence_$(simulation_params["name"])_streamrate.pdf", "w") do file
            writemime(file, "application/pdf", fig)
        end
    end
    PyPlot.close(fig)
end

function plot_SNR(results, simulation_params, precoding_settings, plot_params;
    plot_confidence_interval = false, z_alpha_half = 1.96)

    ### POSTPROCESSING ###
    results_sumrate_mean = [ string(method) => Dict{ASCIIString, Array{Float64}}() for method in simulation_params["precoding_methods"] ]
    results_sumrate_var  = [ string(method) => Dict{ASCIIString, Array{Float64}}() for method in simulation_params["precoding_methods"] ]
    results_minrate_mean = [ string(method) => Dict{ASCIIString, Array{Float64}}() for method in simulation_params["precoding_methods"] ]
    results_minrate_var  = [ string(method) => Dict{ASCIIString, Array{Float64}}() for method in simulation_params["precoding_methods"] ]

    for method_name in intersect(simulation_params["precoding_methods"], keys(plot_params["precoding_methods"]))
        for (result_param, plot_key, plot_legend) in plot_params["precoding_methods"][method_name]
            if isa(result_param, ASCIIString)
                result_name = result_param
                result = results[method_name][result_name]
            else
                (calculator, calculate_from) = result_param
                result_name = string(string(calculator), "_", calculate_from)
                result = calculator(results[method_name][calculate_from])
            end

            # mean: average over drops and sims
            results_sumrate_mean[method_name][result_name] = squeeze(mean(sum(result, 4:5), 1:2), [1,2,4,5])
            results_minrate_mean[method_name][result_name] = squeeze(mean(minimum(sum(result, 5), 4), 1:2), [1,2,4,5])

            # var: average over sims, estimate var over drops
            results_sumrate_var[method_name][result_name] = squeeze(var(mean(sum(result, 4:5), 2), 1), [1,2,4,5])
            results_minrate_var[method_name][result_name] = squeeze(var(mean(minimum(sum(result, 5), 4), 2), 1), [1,2,4,5])
        end
    end

    ### SUM RATE EVOLUTION ###
    fig, ax = plot_methods(simulation_params["Ps_dBm"], results_sumrate_mean,
                results_sumrate_var, simulation_params, plot_params,
                plot_confidence_interval, z_alpha_half)

    #ax[:legend](loc="upper left")
    ax[:set_xlabel]("Transmit power [dBm]")
    ax[:set_ylabel]("Sum rate [bits/s/Hz]")

    fig[:suptitle]("Number of iterations: $(precoding_settings["stop_crit"])")

    if displayable("application/pdf")
        display(fig)
    else
        open("SNR_$(simulation_params["name"])_sumrate.pdf", "w") do file
            writemime(file, "application/pdf", fig)
        end
    end
    PyPlot.close(fig)

    ### MIN RATE EVOLUTION ###
    fig, ax = plot_methods(simulation_params["Ps_dBm"], results_minrate_mean,
                results_minrate_var, simulation_params, plot_params,
                plot_confidence_interval, z_alpha_half)

    #ax[:legend](loc="upper left")
    ax[:set_xlabel]("Transmit power [dBm]")
    ax[:set_ylabel]("Min rate [bits/s/Hz]")

    fig[:suptitle]("Number of iterations: $(precoding_settings["stop_crit"])")

    if displayable("application/pdf")
        display(fig)
    else
        open("SNR_$(simulation_params["name"])_minrate.pdf", "w") do file
            writemime(file, "application/pdf", fig)
        end
    end
    PyPlot.close(fig)
end

function plot_methods(xvals, results_mean, results_var, simulation_params,
    plot_params, plot_confidence_interval, z_alpha_half)

    fig = PyPlot.figure(figsize=plot_params["figsize"])
    ax = fig[:add_subplot](1, 1, 1)

    for method_name in intersect(simulation_params["precoding_methods"], keys(plot_params["precoding_methods"]))
        for (result_param, plot_key, plot_legend) in plot_params["precoding_methods"][method_name]
            if isa(result_param, ASCIIString)
                result_name = result_param
            else
                (calculator, calculate_from) = result_param
                result_name = string(string(calculator), "_", calculate_from)
            end

            ax[:plot](xvals, results_mean[method_name][result_name], plot_key, label=plot_legend)

            if plot_confidence_interval
                facecolor = match(r"[a-z]", plot_key).match
                ax[:fill_between](xvals, results_mean[method_name][result_name] + z_alpha_half*results_var[method_name][result_name]/sqrt(simulation_params["Ndrops"]), results_mean[method_name][result_name] - z_alpha_half*results_var[method_name][result_name]/sqrt(simulation_params["Ndrops"]), facecolor=facecolor, alpha=0.5)
            end
        end
    end

    return (fig, ax)
end
