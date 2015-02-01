function simulate_performance(network::Network,
    simulation_params::SimulationParams,
    precoding_params::PrecodingParams)

    println("-- performance test on $network.")
    Lumberjack.info("Starting performance test.",
        { :network => network, :simulation_params => simulation_params,
          :precoding_params => precoding_params })

    cell_assignment = assign_cells_by_id(network)

    draw_user_drop!(network)
    channel = draw_channel(network)

    # Make sure things are JITed
    for method in simulation_params["precoding_methods"]
        method(channel, network, cell_assignment, precoding_params)
    end
    println("--- JITing @time macro")
    @time 1

    # Run performance test
    for method in simulation_params["precoding_methods"]
        println("--- Testing performance of ", string(method))
        @time for i = 1:simulation_params["Ntest"]; method(channel, network, cell_assignment, precoding_params); end
    end
end
