##########################################################################
# Channels
abstract Channel

immutable SinglecarrierChannel <: Channel
    H::Array{Matrix{Complex128},2}
    Ns::Vector{Int}
    Ms::Vector{Int}
    K::Int
    I::Int
end

immutable MulticarrierChannel <: Channel
    H::Array{Array{Complex128,3},2}
    Ns::Vector{Int}
    Ms::Vector{Int}
    K::Int
    I::Int
    Lc::Int
end

get_channel_gains_dB(channel::SinglecarrierChannel) = 
    Float64[ 20*log10(vecnorm(channel.H[k,i])) for k = 1:channel.K, i = 1:channel.I ]

function get_average_channel_gains_dB(channel::MulticarrierChannel)
    average_channel_gains = zeros(Float64, channel.K, channel.I)

    for k = 1:channel.K
        for i = 1:channel.I
            for l = 1:channel.Lc
                average_channel_gains[k,i] += (1/channel.Lc)*vecnorm(channel.H[k,i,l])^2
            end
        end
    end

    return 10*log10(average_channel_gains)
end
