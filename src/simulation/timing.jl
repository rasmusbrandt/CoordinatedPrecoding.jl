##########################################################################
# Timing for list of methods. loop_over works like in simulate.
function timing(network, simulation_params; loop_over::Symbol=:precoding_methods)
    println("-- timing on $network.")
    Lumberjack.info("Starting timing.",
        [ :network => network, :simulation_params => simulation_params ])

    # Get assignment functions if needed
    _, assignment_method = get_other_method(simulation_params, loop_over)

    # Set initial aux params
    set_initial_aux_params!(simulation_params, network)

    # Ensure that we are not storing all intermediate iterations.
    set_aux_precoding_param!(network, :final_iteration, "output_protocol")

    # Get channel realization to work on
    draw_user_drop!(network)
    channel = draw_channel(network)

    println("--- JITing @time macro")
    @time 1

    if loop_over == :precoding_methods
        # Perform assignment
        assignment_method(channel, network)

        # Make sure things are JITed
        for method in simulation_params["precoding_methods"]
            method(channel, network)
        end

        # Run performance test
        for method in simulation_params["precoding_methods"]
            println("--- Testing performance of ", string(method))
            @time for i = 1:simulation_params["Ntest"]; method(channel, network); end
        end
    elseif loop_over == :assignment_methods
        # Make sure things are JITed
        for method in simulation_params["assignment_methods"]
            method(channel, network)
        end

        # Run performance test
        for method in simulation_params["assignment_methods"]
            println("--- Testing performance of ", string(method))
            @time for i = 1:simulation_params["Ntest"]; method(channel, network); end
        end
    end
end
