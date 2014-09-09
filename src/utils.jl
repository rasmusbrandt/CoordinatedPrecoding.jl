# Some of the things here might be weird. To be improved, or removed!


# For some reason, the base library does not seem to handle Hermitian matrices perfectly yet. Pull request?
+(A::Hermitian{Complex{Float64}}, B::Hermitian{Complex{Float64}}) = Hermitian(A.S + B.S)


# need to fix these!!!
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
