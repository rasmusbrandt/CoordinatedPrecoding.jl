using CoordinatedPrecoding
using Base.Test

function build_identity_channel(N, K)
    coefs = Array(Matrix{Complex128}, K, K)

    for k = 1:K
        for i = 1:K
            coefs[k,i] = (1/sqrt(2))*(eye(N, N) + im*eye(N, N))
        end
    end

    SinglecarrierChannel(coefs, N*ones(Int, K), N*ones(Int, K), K, K)
end

N = 2; K = 3
@test_approx_eq 10*log10(N*ones(Float64, K, K)) get_channel_gains_dB(build_identity_channel(N, K))
