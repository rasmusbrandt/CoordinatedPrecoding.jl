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

# All of these things should ideally be in Base. Investigate and pull request?
+(A::Hermitian{Complex{Float64}}, B::Hermitian{Complex{Float64}}) = Hermitian(A.S + B.S)

function +(A::Hermitian{Complex{Float64}}, B::Array{Float64,2})
    if ishermitian(B)
        return Hermitian(A.S + B)
    else
        return A.S + B
    end
end
+(B::Array{Float64,2}, A::Hermitian{Complex{Float64}}) = +(A::Hermitian{Complex{Float64}}, B::Array{Float64,2})

function -(A::Hermitian{Complex{Float64}}, B::Array{Complex{Float64},2})
    if ishermitian(B)
        return Hermitian(A.S - B)
    else
        return A.S - B
    end
end
-(B::Array{Complex{Float64},2}, A::Hermitian{Complex{Float64}}) = -(A::Hermitian{Complex{Float64}}, B::Array{Complex{Float64},2})

.*(a::Float64, B::Hermitian{Complex{Float64}}) = a.*B.S

logdet(A::Hermitian{Complex{Float64}}) = Base.logdet(A.S)
