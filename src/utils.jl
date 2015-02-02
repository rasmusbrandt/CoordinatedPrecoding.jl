##########################################################################
# General utilities
function clean_simulation_params_for_jld(simulation_params)
    cleaned_precoding_methods = Array(ASCIIString, 0)
    for method in simulation_params["precoding_methods"]
        push!(cleaned_precoding_methods, string(method))
    end

    cleaned_aux_precoding_params = Dict{ASCIIString, Any}()
    for (key, val) in simulation_params["aux_precoding_params"]
        if isa(val, Function)
            cleaned_aux_precoding_params[key] = string(val)
        else
            cleaned_aux_precoding_params[key] = val
        end
    end

    cleaned_simulation_params = copy(simulation_params)
    cleaned_simulation_params["precoding_methods"] = cleaned_precoding_methods
    cleaned_simulation_params["aux_precoding_params"] = cleaned_aux_precoding_params

    return cleaned_simulation_params
end

##########################################################################
# Hermitian utility functions. These should probably be in Base. Contribute?
+(A::Hermitian{Complex128}, B::Hermitian{Complex128}) = Hermitian(A.S + B.S)
+(B::Matrix{Float64}, A::Hermitian{Complex128}) = +(A, B)
+(A::Hermitian{Complex128}, B::Matrix{Float64}) = A.S + B
+(A::Hermitian{Complex128}, B::Matrix{Complex128}) = A.S + B

-(A::Hermitian{Complex128}, B::Matrix{Complex128}) = A.S - B
-(B::Array{Complex128,2}, A::Hermitian{Complex128}) = -(A, B)
-(A::Hermitian{Complex128}, B::Matrix{Float64}) = A.S - B

.*(a::Float64, B::Hermitian{Complex128}) = Hermitian(a.*B.S)

import Base.logdet, Base.diag
logdet(A::Hermitian{Complex128}) = logdet(A.S)
diag(A::Hermitian{Complex128}) = diag(A.S)
