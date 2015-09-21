##########################################################################
# General utilities
stringify_method(f::Function) = string(f.env.name)
@compat stringify_method(f::AbstractString) = f

function clean_simulation_params_for_jld(simulation_params)
    cleaned_simulation_params = copy(simulation_params)

    # Precoding method function pointers
    if haskey(simulation_params, "precoding_methods")
        cleaned_precoding_methods = Array(ASCIIString, 0)
        for method in simulation_params["precoding_methods"]
            push!(cleaned_precoding_methods, stringify_method(method))
        end
        cleaned_simulation_params["precoding_methods"] = cleaned_precoding_methods
    end
    if haskey(simulation_params, "precoding_methods_nosweep")
        cleaned_precoding_methods_nosweep = Array(ASCIIString, 0)
        for method in simulation_params["precoding_methods_nosweep"]
            push!(cleaned_precoding_methods_nosweep, stringify_method(method))
        end
        cleaned_simulation_params["precoding_methods_nosweep"] = cleaned_precoding_methods_nosweep
    end

    # Assignment method function pointers
    if haskey(simulation_params, "assignment_methods")
        cleaned_assignment_methods = Array(ASCIIString, 0)
        for method in simulation_params["assignment_methods"]
            push!(cleaned_assignment_methods, stringify_method(method))
        end
        cleaned_simulation_params["assignment_methods"] = cleaned_assignment_methods
    end

    # Independent variable function pointer
    if haskey(simulation_params, "independent_variable")
        cleaned_independent_variable = (string(simulation_params["independent_variable"][1]), simulation_params["independent_variable"][2])
        cleaned_simulation_params["independent_variable"] = cleaned_independent_variable
    end

    # Auxiliary independent variable function pointers
    if haskey(simulation_params, "aux_independent_variables")
        cleaned_aux_independent_variables = Array((@compat Tuple{ASCIIString, Any}), 0)
        for aux_idp in simulation_params["aux_independent_variables"]
            push!(cleaned_aux_independent_variables, (string(aux_idp[1]), aux_idp[2]))
        end
        cleaned_simulation_params["aux_independent_variables"] = cleaned_aux_independent_variables
    end

    return cleaned_simulation_params
end

##########################################################################
# Parameter dict defaultization
macro defaultize_param!(params, key, default_val)
    # Note that we are intentionally violating macro hygiene here...
    return esc(:(haskey($params, $key) || setindex!($params, $default_val, $key)))
end
