##########################################################################
# Timing for list of methods. loop_over works like in simulate.
function timing(network, simulation_params; loop_over::Symbol=:precoding_methods)
    println("-- timing on $network.")
    Lumberjack.info("Starting timing.", @compat Dict(:network => network, :simulation_params => simulation_params))

    # Get assignment method if needed
    assignment_method() = nothing
    if loop_over == :precoding_methods
        if haskey(simulation_params, "assignment_methods")
            if length(simulation_params["assignment_methods"]) > 1
                Lumberjack.warn("Looping over precoding methods: will only use first assignment method provided.")
            end
            assignment_method(channel, network) = simulation_params["assignment_methods"][1](channel, network)
        else
            assignment_method(channel, network) = IDCellAssignment!(channel, network)
        end
    end

    # Set initial aux params
    set_initial_aux_params!(simulation_params, network)

    # Ensure that we are not storing all intermediate iterations.
    set_aux_precoding_param!(network, :final_iteration, "output_protocol")

    # Get channel realization to work on
    draw_user_drop!(network)
    channel = draw_channel(network)

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
