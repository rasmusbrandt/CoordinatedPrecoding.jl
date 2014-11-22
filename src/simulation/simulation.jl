function simulate_convergence(network, simulation_params, precoding_settings)
    println("-- simulate_convergence on $network.")
    println("--- Ndrops: $(simulation_params["Ndrops"]), Nsim: $(simulation_params["Nsim"]).")
    Lumberjack.info("Starting convergence simulation.",
        { :network => network, :simulation_params => simulation_params,
          :precoding_settings => precoding_settings })

    cell_assignment = assign_cells_by_id(network)

    # Ensure that we are using output protocol 1, because we want to store all
    # intermediate iterations.
    precoding_settings["output_protocol"] = 1

    results = [ string(method) => Dict{ASCIIString, Array{Float64}}() for method in simulation_params["precoding_methods"] ]
    result_ranges = [ string(method) => Dict{ASCIIString, Array{UnitRange{Int64},1}}() for method in simulation_params["precoding_methods"] ]

    tic()
    for Ndrops_ind = 1:simulation_params["Ndrops"]
        if Ndrops_ind == simulation_params["Ndrops"]
            println(Ndrops_ind, ".")
        else
            print(Ndrops_ind, ", ")
        end

        draw_user_drop!(network)

        for Nsim_ind = 1:simulation_params["Nsim"]
            channel = draw_channel(network)

            for method in simulation_params["precoding_methods"]
                method_results = method(channel, network, cell_assignment, precoding_settings)

                for (result_name, result) in method_results
                    if !haskey(results[string(method)], result_name)
                        # Create storage container and remember dimensions
                        results[string(method)][result_name] = Array(Float64, simulation_params["Ndrops"], simulation_params["Nsim"], size(result)...)
                        result_ranges[string(method)][result_name] = [ 1:s for s in size(result) ]
                    end

                    results[string(method)][result_name][Ndrops_ind, Nsim_ind, result_ranges[string(method)][result_name]...] = result
                end
            end
        end
    end
    t = toq(); println("--- elapsed time: ", t/60, " minutes")

    return results
end

function simulate_SNR(network, simulation_params, precoding_settings)
    println("-- simulate_SNR on $network.")
    println("--- Ndrops: $(simulation_params["Ndrops"]), Nsim: $(simulation_params["Nsim"]).")
    Lumberjack.info("Starting SNR simulation.",
        { :network => network, :simulation_params => simulation_params,
          :precoding_settings => precoding_settings })

    cell_assignment = assign_cells_by_id(network)

    # Ensure that we are using output protocol 2, so we don't have to store all
    # the intermediate iterations.
    precoding_settings["output_protocol"] = 2

    results = [ string(method) => Dict{ASCIIString, Array{Float64}}() for method in simulation_params["precoding_methods"] ]
    result_ranges = [ string(method) => Dict{ASCIIString, Array{UnitRange{Int64},1}}() for method in simulation_params["precoding_methods"] ]

    tic()
    for Ndrops_ind = 1:simulation_params["Ndrops"]
        if Ndrops_ind == simulation_params["Ndrops"]
            println(Ndrops_ind, ".")
        else
            print(Ndrops_ind, ", ")
        end

        draw_user_drop!(network)

        for Nsim_ind = 1:simulation_params["Nsim"]
            channel = draw_channel(network)

            for P_dBm_ind in 1:length(simulation_params["Ps_dBm"])
                set_transmit_powers!(network, 10^(simulation_params["Ps_dBm"][P_dBm_ind]/10))

                for method in simulation_params["precoding_methods"]
                    method_results = method(channel, network, cell_assignment, precoding_settings)

                    for (result_name, result) in method_results
                        if !haskey(results[string(method)], result_name)
                            # Create storage container and remember dimensions
                            results[string(method)][result_name] = Array(Float64, simulation_params["Ndrops"], simulation_params["Nsim"], length(simulation_params["Ps_dBm"]), size(result)...)
                            result_ranges[string(method)][result_name] = [ 1:s for s in size(result) ]
                        end

                        results[string(method)][result_name][Ndrops_ind, Nsim_ind, P_dBm_ind, result_ranges[string(method)][result_name]...] = result
                    end
                end
            end
        end
    end
    t = toq(); println("--- elapsed time: ", t/60, " minutes")

    return results
end

function perform_performancetest(network, simulation_params, precoding_settings)
    println("-- performance test on $network.")
    Lumberjack.info("Starting performance test.",
        { :network => network, :simulation_params => simulation_params,
          :precoding_settings => precoding_settings })

    cell_assignment = assign_cells_by_id(network)

    draw_user_drop!(network)
    channel = draw_channel(network)

    # Make sure things are JITed
    for method in simulation_params["precoding_methods"]
        method(channel, network, cell_assignment, precoding_settings)
    end
    println("--- JITing @time macro")
    @time 1

    # Run performance test
    for method in simulation_params["precoding_methods"]
        println("--- Testing performance of ", string(method))
        @time for i = 1:simulation_params["Ntest"]; method(channel, network, cell_assignment, precoding_settings); end
    end
end
