##########################################################################
# General utilities
function clean_simulation_params_for_jld(simulation_params)
    cleaned_simulation_params = copy(simulation_params)

    # Precoding method function pointers
    if haskey(simulation_params, "precoding_methods")
        cleaned_precoding_methods = Array(ASCIIString, 0)
        for method in simulation_params["precoding_methods"]
            push!(cleaned_precoding_methods, string(method))
        end
        cleaned_simulation_params["precoding_methods"] = cleaned_precoding_methods
    end

    # Precoding params function pointers
    if haskey(simulation_params, "aux_precoding_params")
        cleaned_aux_precoding_params = Dict{ASCIIString, Any}()
        for (key, val) in simulation_params["aux_precoding_params"]
            if isa(val, Function)
                cleaned_aux_precoding_params[key] = string(val)
            else
                cleaned_aux_precoding_params[key] = val
            end
        end
        cleaned_simulation_params["aux_precoding_params"] = cleaned_aux_precoding_params
    end

    # Independent variable function pointer
    if haskey(simulation_params, "independent_variable")
        cleaned_independent_variable = (string(simulation_params["independent_variable"][1]), simulation_params["independent_variable"][2])
        cleaned_simulation_params["independent_variable"] = cleaned_independent_variable
    end

    # Auxiliary independent variable function pointers
    if haskey(simulation_params, "aux_independent_variables")
        cleaned_aux_independent_variables = Array((ASCIIString, Any), 0)
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

##########################################################################
# Hermitian utility functions. These should probably be in Base. Contribute?
+(A::Hermitian{Complex128}, B::Hermitian{Complex128}) = Hermitian(A.S + B.S)
+(B::Matrix{Float64}, A::Hermitian{Complex128}) = +(A, B)
+(A::Hermitian{Complex128}, B::Matrix{Float64}) = A.S + B
+(A::Hermitian{Complex128}, B::Matrix{Complex128}) = A.S + B

-(A::Hermitian{Complex128}, B::Matrix{Complex128}) = A.S - B
-(B::Array{Complex128, 2}, A::Hermitian{Complex128}) = -(A, B)
-(A::Hermitian{Complex128}, B::Matrix{Float64}) = A.S - B

.*(a::Float64, B::Hermitian{Complex128}) = Hermitian(a.*B.S)

import Base.logdet, Base.diag
logdet(A::Hermitian{Complex128}) = logdet(A.S)
diag(A::Hermitian{Complex128}) = diag(A.S)
