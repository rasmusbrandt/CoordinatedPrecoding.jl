using CoordinatedPrecoding
using PyCall; @pyimport matplotlib.pyplot as plt

##########################################################################
# Simulation and visualization functions
function simulate_SNR(simulation_params, precoding_settings)
    network = setup_triangular3site_network(simulation_params["Kc"],
        simulation_params["N"], simulation_params["M"],
        no_streams=simulation_params["d"])

    results1 = Array(Float64, length(simulation_params["Ps_dBm"]), get_no_MSs(network), simulation_params["d"], simulation_params["Ndrops"], simulation_params["Nshadow"], simulation_params["Nsim"])
    results2 = Array(Float64, length(simulation_params["Ps_dBm"]), get_no_MSs(network), simulation_params["d"], simulation_params["Ndrops"], simulation_params["Nshadow"], simulation_params["Nsim"])
    results3 = Array(Float64, length(simulation_params["Ps_dBm"]), get_no_MSs(network), simulation_params["d"], simulation_params["Ndrops"], simulation_params["Nshadow"], simulation_params["Nsim"])

    for Ndrops_ind = 1:simulation_params["Ndrops"]
        println("Drop $Ndrops_ind")
        drop_users_randomly!(network)

        for Nshadow_ind = 1:simulation_params["Nshadow"]
            draw_shadow_fading!(network)

            for Nsim_ind = 1:simulation_params["Nsim"]
                H = draw_channel(network)

                for P_dB_ind in 1:length(simulation_params["Ps_dBm"])
                    set_transmit_powers!(network, 10^(simulation_params["Ps_dBm"][P_dB_ind]/10))

                    results1[P_dB_ind, :, :, Ndrops_ind, Nshadow_ind, Nsim_ind] = Shi2011_WMMSE(H, network, precoding_settings)[:,:,end]
                    results2[P_dB_ind, :, :, Ndrops_ind, Nshadow_ind, Nsim_ind] = Gomadam2008_MaxSINR(H, network, precoding_settings)[:,:,end]
                    results3[P_dB_ind, :, :, Ndrops_ind, Nshadow_ind, Nsim_ind] = Komulainen2013_WMMSE(H, network, precoding_settings)[:,:,end]
                end
            end
        end
    end

    [ "Shi2011_WMMSE" => results1, 
      "Gomadam2008_MaxSINR" => results2,
      "Komulainen2013_WMMSE" => results3,
    ]
end

function plot_SNR(results, simulation_params, precoding_settings)
    fig = plt.figure(figsize=(8,4))
    ax = fig[:add_subplot](1,1,1)

    ax[:plot](simulation_params["Ps_dBm"], squeeze(mean(sum(results["Shi2011_WMMSE"], [2,3]), [4,5,6]), [2,3,4,5,6]), "b-", label="Shi2011_WMMSE")
    ax[:plot](simulation_params["Ps_dBm"], squeeze(mean(sum(results["Gomadam2008_MaxSINR"], [2,3]), [4,5,6]), [2,3,4,5,6]), "r-", label="Gomadam2008_MaxSINR")
    ax[:plot](simulation_params["Ps_dBm"], squeeze(mean(sum(results["Komulainen2013_WMMSE"], [2,3]), [4,5,6]), [2,3,4,5,6]), "b-.", label="Komulainen2013_WMMSE")

    ax[:legend](loc="best")
    ax[:set_xlabel]("Transmit power [dBm]")
    ax[:set_ylabel]("Sum rate rate [bits/s/Hz]")

    fig[:suptitle]("Number of iterations: $(precoding_settings["stop_crit"])")

    fig[:savefig]("triangular3site_example.png", dpi=125)
end

##########################################################################
# Settings
srand(211438269237)

const simulations_params = {
    "Kc" => 2, "N" => 3, "M" => 3,
    "Ndrops" => 2, "Nshadow" => 2, "Nsim" => 2,
    "Ps_dBm" => -10:3:20,
    "d" => 2
}
const precoding_settings = { "stop_crit" => 50 }

results = simulate_SNR(simulations_params, precoding_settings)
plot_SNR(results, simulations_params, precoding_settings)
