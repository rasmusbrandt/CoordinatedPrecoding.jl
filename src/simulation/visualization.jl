##########################################################################
# Postprocess results from the simulate function.
function postprocess(raw_results, simulation_params, plot_params)
    methods = get_methods_to_plot(simulation_params, plot_params)

    # Main independent variable
    idp_vals_length = length(simulation_params["independent_variable"][2])

    # Auxiliary independent variables
    if haskey(simulation_params, "aux_independent_variables")
        Naux = length(simulation_params["aux_independent_variables"])
        aux_idp_vals_length = length(simulation_params["aux_independent_variables"][1][2])
        for n = 2:Naux
            aux_idp_vals_length == length(simulation_params["aux_independent_variables"][n][2]) ? nothing : Lumberjack.error("Auxiliary independent variable vectors must have equal length.")
        end
    else
        Naux = 0; aux_idp_vals_length = 1
    end

    # Compact raw results into result matrices
    Ndrops = size(raw_results.simulation_results, 1); Nsim = size(raw_results.simulation_results, 2)
    results = [ method => Dict{ASCIIString, Array{Float64}}() for method in methods ]
    for method_name in methods
        for (result_param,) in plot_params["methods"][method_name]
            if isa(result_param, ASCIIString)
                result_name = result_param

                result_dimensions = size(raw_results[1, 1, 1, 1][method_name][result_name])
            else
                (calculator, calculate_from) = result_param
                result_name = string(string(calculator), "_", calculate_from)

                result_dimensions = size(raw_results[1, 1, 1, 1][method_name][calculate_from])
            end
            result_ranges = [ 1:s for s in result_dimensions ]
            results[method_name][result_name] = Array(Float64, Ndrops, Nsim, idp_vals_length, aux_idp_vals_length, result_dimensions...)
            for Ndrops_idx = 1:Ndrops; for Nsim_idx = 1:Nsim; for idp_vals_idx = 1:idp_vals_length; for aux_idp_vals_idx = 1:aux_idp_vals_length
                if isa(result_param, ASCIIString)
                    result = raw_results[Ndrops_idx, Nsim_idx, idp_vals_idx, aux_idp_vals_idx][method_name][result_name]
                else
                    result = calculator(raw_results[Ndrops_idx, Nsim_idx, idp_vals_idx, aux_idp_vals_idx][method_name][calculate_from], idp_vals_idx)
                end
                results[method_name][result_name][Ndrops_idx, Nsim_idx, idp_vals_idx, aux_idp_vals_idx, result_ranges...] = result
            end; end; end; end
        end
    end

    # Calculate statistics on result matrices
    plot_params["objective"] == :sum && (objective = (r -> sum(r, 5:6)))
    plot_params["objective"] == :mean && (objective = (r -> mean(sum(r, 6), 5)))
    plot_params["objective"] == :min && (objective = (r -> minimum(sum(r, 6), 5)))
    plot_params["objective"] == :none && (objective = (r -> r))

    results_mean = [ method => Dict{ASCIIString, Array{Float64}}() for method in methods ]
    results_var = [ method => Dict{ASCIIString, Array{Float64}}() for method in methods ]
    for method_name in methods
        for (result_param,) in plot_params["methods"][method_name]
            if isa(result_param, ASCIIString)
                result_name = result_param
            else
                (calculator, calculate_from) = result_param
                result_name = string(string(calculator), "_", calculate_from)
            end

            # mean: average over drops and sims
            results_mean[method_name][result_name] = squeeze(mean(objective(results[method_name][result_name]), 1:2), (1,2,5,6))

            # var: average over sims, estimate var over drops
            results_var[method_name][result_name] = squeeze(var(mean(objective(results[method_name][result_name]), 2), 1), (1,2,5,6))
        end
    end

    return results, results_mean, results_var
end

##########################################################################
# Plot results from the postprocess function.
function plot(processed_results, simulation_params, plot_params)
    results_mean = processed_results[2]; results_var = processed_results[3]

    fig = PyPlot.figure(;plot_params["figure"]...)
    ax = fig[:add_subplot](1, 1, 1; plot_params["axes"]...)

    plot_methods(ax, results_mean, results_var, simulation_params, plot_params;
        xvals=(haskey(plot_params, "xvals") ? plot_params["xvals"] : simulation_params["independent_variable"][2]))
    haskey(plot_params, "legend") && ax[:legend](;plot_params["legend"]...)

    if displayable("application/pdf") || displayable("image/png")
        display(fig)
    else
        open("$(simulation_params["simulation_name"])_$(plot_params["plot_name"]).pdf", "w") do file
            writemime(file, "application/pdf", fig)
        end
    end
    PyPlot.close(fig)
end

##########################################################################
# Post-process results from the simulate_precoding_convergence function.
function postprocess_precoding_convergence(raw_results, simulation_params, plot_params)
    methods = get_methods_to_plot(simulation_params, plot_params)

    # Auxiliary independent variables
    if haskey(simulation_params, "aux_independent_variables")
        Naux = length(simulation_params["aux_independent_variables"])
        aux_idp_vals_length = length(simulation_params["aux_independent_variables"][1][2])
        for n = 2:Naux
            aux_idp_vals_length == length(simulation_params["aux_independent_variables"][n][2]) ? nothing : Lumberjack.error("Auxiliary independent variable vectors must have equal length.")
        end
    else
        Naux = 0; aux_idp_vals_length = 1
    end

    # Compact raw results into result matrices
    Ndrops = size(raw_results.simulation_results, 1); Nsim = size(raw_results.simulation_results, 2)
    results = [ method => Dict{ASCIIString, Array{Float64}}() for method in methods ]
    for method_name in methods
        for (result_param,) in plot_params["methods"][method_name]
            if isa(result_param, ASCIIString)
                result_name = result_param

                result_dimensions = size(raw_results[1, 1, 1][method_name][result_name])
            else
                (calculator, calculate_from) = result_param
                result_name = string(string(calculator), "_", calculate_from)

                result_dimensions = size(raw_results[1, 1, 1][method_name][calculate_from])
            end

            result_ranges = [ 1:s for s in result_dimensions ]
            results[method_name][result_name] = Array(Float64, Ndrops, Nsim, aux_idp_vals_length, result_dimensions...)
            for Ndrops_idx = 1:Ndrops; for Nsim_idx = 1:Nsim; for aux_idp_vals_idx = 1:aux_idp_vals_length
                if isa(result_param, ASCIIString)
                    result = raw_results[Ndrops_idx, Nsim_idx, aux_idp_vals_idx][method_name][result_name]
                else
                    result = calculator(raw_results[Ndrops_idx, Nsim_idx, aux_idp_vals_idx][method_name][calculate_from])
                end
                results[method_name][result_name][Ndrops_idx, Nsim_idx, aux_idp_vals_idx, result_ranges...] = result
            end; end; end
        end
    end

    # Calculate statistics on result matrices
    plot_params["objective"] == :sum && (objective = (r -> sum(r, 4:5)))
    plot_params["objective"] == :mean && (objective = (r -> mean(sum(r, 5), 4)))
    plot_params["objective"] == :min && (objective = (r -> minimum(sum(r, 5), 4)))
    plot_params["objective"] == :none && (objective = (r -> r))

    results_mean = [ method => Dict{ASCIIString, Array{Float64}}() for method in methods ]
    results_var = [ method => Dict{ASCIIString, Array{Float64}}() for method in methods ]
    for method_name in methods
        for (result_param,) in plot_params["methods"][method_name]
            if isa(result_param, ASCIIString)
                result_name = result_param
            else
                (calculator, calculate_from) = result_param
                result_name = string(string(calculator), "_", calculate_from)
            end

            # mean: average over drops and sims
            results_mean[method_name][result_name] = transpose(squeeze(mean(objective(results[method_name][result_name]), 1:2), (1,2,4,5)))

            # var: average over sims, estimate var over drops
            results_var[method_name][result_name] = transpose(squeeze(var(mean(objective(results[method_name][result_name]), 2), 1), (1,2,4,5)))
        end
    end

    return results, results_mean, results_var
end

##########################################################################
# Post-process results from the simulate_assignment_convergence function.
function postprocess_assignment_convergence(raw_results, simulation_params, plot_params)
    methods = get_methods_to_plot(simulation_params, plot_params)

    # Auxiliary independent variables
    if haskey(simulation_params, "aux_independent_variables")
        Naux = length(simulation_params["aux_independent_variables"])
        aux_idp_vals_length = length(simulation_params["aux_independent_variables"][1][2])
        for n = 2:Naux
            aux_idp_vals_length == length(simulation_params["aux_independent_variables"][n][2]) ? nothing : Lumberjack.error("Auxiliary independent variable vectors must have equal length.")
        end
    else
        Naux = 0; aux_idp_vals_length = 1
    end

    # Compact raw results into result matrices
    Ndrops = size(raw_results.simulation_results, 1)
    results = [ method => Dict{ASCIIString, Array{Float64}}() for method in methods ]
    for method_name in methods
        for (result_param,) in plot_params["methods"][method_name]
            if isa(result_param, ASCIIString)
                result_name = result_param

                result_dimensions = size(raw_results[1, 1, 1][method_name][result_name])
            else
                (calculator, calculate_from) = result_param
                result_name = string(string(calculator), "_", calculate_from)

                result_dimensions = size(raw_results[1, 1, 1][method_name][calculate_from])
            end

            result_ranges = [ 1:s for s in result_dimensions ]
            results[method_name][result_name] = Array(Float64, Ndrops, 1, aux_idp_vals_length, result_dimensions...)
            for Ndrops_idx = 1:Ndrops; for aux_idp_vals_idx = 1:aux_idp_vals_length
                if isa(result_param, ASCIIString)
                    result = raw_results[Ndrops_idx, 1, aux_idp_vals_idx][method_name][result_name]
                else
                    result = calculator(raw_results[Ndrops_idx, 1, aux_idp_vals_idx][method_name][calculate_from])
                end
                results[method_name][result_name][Ndrops_idx, 1, aux_idp_vals_idx, result_ranges...] = result
            end; end
        end
    end

    # Calculate statistics on result matrices
    plot_params["objective"] == :sum && (objective = (r -> sum(r, 4:5)))
    plot_params["objective"] == :mean && (objective = (r -> mean(sum(r, 5), 4)))
    plot_params["objective"] == :min && (objective = (r -> minimum(sum(r, 5), 4)))
    plot_params["objective"] == :none && (objective = (r -> r))

    results_mean = [ method => Dict{ASCIIString, Array{Float64}}() for method in methods ]
    results_var = [ method => Dict{ASCIIString, Array{Float64}}() for method in methods ]
    for method_name in methods
        for (result_param,) in plot_params["methods"][method_name]
            if isa(result_param, ASCIIString)
                result_name = result_param
            else
                (calculator, calculate_from) = result_param
                result_name = string(string(calculator), "_", calculate_from)
            end

            # mean: average over drops and sims
            results_mean[method_name][result_name] = transpose(squeeze(mean(objective(results[method_name][result_name]), 1:2), (1,2,4,5)))

            # var: average over sims, estimate var over drops
            results_var[method_name][result_name] = transpose(squeeze(var(mean(objective(results[method_name][result_name]), 2), 1), (1,2,4,5)))
        end
    end

    return results, results_mean, results_var
end

##########################################################################
# Plots results from the postprocess_precoding_convergence function.
function plot_precoding_convergence(processed_results, simulation_params, plot_params; user_plots::Bool=true)
    methods = get_methods_to_plot(simulation_params, plot_params)

    results = processed_results[1]
    results_mean = processed_results[2]
    results_var = processed_results[3]

    ### SYSTEM-LEVEL OBJECTIVE ###
    fig = PyPlot.figure(;plot_params["figure"]...)
    ax = fig[:add_subplot](1, 1, 1; plot_params["axes"]...)

    if haskey(plot_params, "xvals")
        plot_methods(ax, results_mean, results_var, simulation_params, plot_params;
            xvals=plot_params["xvals"])
    else
        plot_methods(ax, results_mean, results_var, simulation_params, plot_params)
    end
    haskey(plot_params, "legend") && ax[:legend](;plot_params["legend"]...)

    if displayable("application/pdf") || displayable("image/png")
        display(fig)
    else
        open("$(simulation_params["simulation_name"])_$(plot_params["plot_name"]).pdf", "w") do file
            writemime(file, "application/pdf", fig)
        end
    end
    PyPlot.close(fig)

    if user_plots
        ### USER UTILITIES ###
        K_method_name = collect(keys(results))[1]
        K_result_name = collect(keys(results[K_method_name]))[1]
        K = size(results[K_method_name][K_result_name], 4) # ugly shit. This is really not stable.

        fig = PyPlot.figure(figsize=(6*K,3*length(methods)))
        subplot_ind = 1

        for method_name in methods
            for k = 1:K
                ax = fig[:add_subplot](length(methods), K, subplot_ind); subplot_ind += 1

                for (result_param, line_params) in plot_params["methods"][method_name]
                    if isa(result_param, ASCIIString)
                        result = results[method_name][result_param]
                    else
                        (calculator, calculate_from) = result_param
                        result = calculator(results[method_name][calculate_from])
                    end

                    ax[:plot](1:size(result, 6),
                              transpose(squeeze(mean(sum(result[:,:,:,k,:,:], 5), 1:2), (1,2,4,5)));
                              line_params...)
                end

                if k == 1
                    ax[:set_ylabel](method_name)
                end
                if subplot_ind-1 <= K
                    ax[:set_title]("User $k")
                end
                if subplot_ind-1 > (K-1)*length(methods)
                    ax[:set_xlabel]("Iteration")
                end
            end
        end

        if displayable("application/pdf") || displayable("image/png")
            display(fig)
        else
            open("$(simulation_params["simulation_name"])_$(plot_params["plot_name"])_peruser.pdf", "w") do file
                writemime(file, "application/pdf", fig)
            end
        end
        PyPlot.close(fig)

        ### STREAM UTILITIES ###
        fig = PyPlot.figure(figsize=(6*simulation_params["d"],3*K))
        subplot_ind = 1

        for k = 1:K; for n = 1:simulation_params["d"]
            ax = fig[:add_subplot](K, simulation_params["d"], subplot_ind); subplot_ind += 1

            for method_name in methods
                for (result_param, line_params) in plot_params["methods"][method_name]
                    if isa(result_param, ASCIIString)
                        result = results[method_name][result_param]
                    else
                        (calculator, calculate_from) = result_param
                        result = calculator(results[method_name][calculate_from])
                    end

                    ax[:plot](1:size(result, 6),
                              transpose(squeeze(mean(result[:,:,:,k,n,:], 1:2), (1,2,4,5)));
                              line_params...)
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
        end; end

        if displayable("application/pdf") || displayable("image/png")
            display(fig)
        else
            open("$(simulation_params["simulation_name"])_$(plot_params["plot_name"])_perstream.pdf", "w") do file
                writemime(file, "application/pdf", fig)
            end
        end
        PyPlot.close(fig)
    end
end

##########################################################################
# Plots results from the postprocess_assignment_convergence function.
plot_assignment_convergence(processed_results, simulation_params, plot_params; user_plots::Bool=true) =
    plot_precoding_convergence(processed_results, simulation_params, plot_params, user_plots=user_plots)

##########################################################################
# Helper methods

# Plot lines for all methods into an axes
function plot_methods(ax, results_mean, results_var, simulation_params, plot_params; xvals=[])
    methods = get_methods_to_plot(simulation_params, plot_params)

    for method_name in methods
        for (result_param, line_params) in plot_params["methods"][method_name]
            if isa(result_param, ASCIIString)
                result_name = result_param
            else
                (calculator, calculate_from) = result_param
                result_name = string(string(calculator), "_", calculate_from)
            end

            xvals = (xvals != []) ? xvals : 1:size(results_mean[method_name][result_name], 1)
            ax[:plot](xvals, results_mean[method_name][result_name]; line_params...)

            if haskey(plot_params, "confidence_interval_zalpha2")
                for n = 1:size(results_mean[method_name][result_name], 2)
                    m = results_mean[method_name][result_name][:,n]
                    s = sqrt(results_var[method_name][result_name][:,n])
                    e = plot_params["confidence_interval_zalpha2"]*s/sqrt(simulation_params["Ndrops"])
                    ax[:fill_between](xvals, m + e, m - e; alpha=0.5, color=line_params[:color])
                end
            end
        end
    end
end

# We can only plot the methods that have been simulated (and thus specified
# in simulation_params), and those we have plotting parameters for (and
# thus specified in plot_params).
function get_methods_to_plot(simulation_params, plot_params)
    precoding_methods = get(simulation_params, "precoding_methods", Function[])
    assignment_methods = get(simulation_params, "assignment_methods", Function[])

    return intersect(map(stringify_method, union(precoding_methods, assignment_methods)), keys(plot_params["methods"]))
end
