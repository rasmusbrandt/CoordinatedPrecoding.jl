##########################################################################
# General utilities
function clean_simulation_params_for_jld(simulation_params)
    # Remove all function pointers from precoding_methods
    cleaned_precoding_methods = Array(ASCIIString, 0)
    for method in simulation_params["precoding_methods"]
        push!(cleaned_precoding_methods, string(method))
    end

    cleaned_simulation_params = copy(simulation_params)
    cleaned_simulation_params["precoding_methods"] = cleaned_precoding_methods

    return cleaned_simulation_params
end

function clean_precoding_settings_for_jld(precoding_settings)
    # Remove all function pointers from precoding_methods
    cleaned_precoding_settings = Dict{ASCIIString, Any}()
    for (key, val) in precoding_settings
        if isa(val, Function)
            cleaned_precoding_settings[key] = string(val)
        else
            cleaned_precoding_settings[key] = val
        end
    end

    return cleaned_precoding_settings
end

##########################################################################
# Hermitian utility functions. These should probably be in Base. Contribute?
+(A::Hermitian{Complex128}, B::Hermitian{Complex128}) = Hermitian(A.S + B.S)
+(A::Hermitian{Complex128}, B::Matrix{Float64}) = A.S + B
+(B::Matrix{Float64}, A::Hermitian{Complex128}) = +(A, B)

-(A::Hermitian{Complex128}, B::Array{Complex128,2}) = A.S - B
-(B::Array{Complex128,2}, A::Hermitian{Complex128}) = -(A, B)

.*(a::Float64, B::Hermitian{Complex128}) = Hermitian(a.*B.S)

import Base.logdet, Base.diag
logdet(A::Hermitian{Complex128}) = logdet(A.S)
diag(A::Hermitian{Complex128}) = diag(A.S)
