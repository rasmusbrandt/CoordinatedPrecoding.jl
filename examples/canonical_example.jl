using CoordinatedPrecoding
using PyCall; @pyimport matplotlib.pyplot as plt

##########################################################################
# Simulation and visualization functions
function simulate_single_canonical_realization(simulation_params, precoding_settings)
    network = setup_interference_channel(simulation_params["K"],
        simulation_params["N"], simulation_params["M"],
        transmit_power=10^(simulation_params["P_dB"]/10),
        receiver_noise_power=10^(simulation_params["sigma2_dB"]/10),
        no_streams=simulation_params["d"])

    channel = draw_channel(network)

    user_rates1 = Shi2011_WMMSE(channel, network, precoding_settings)
    user_rates2 = Gomadam2008_MaxSINR(channel, network, precoding_settings)
    user_rates3 = Komulainen2013_WMMSE(channel, network, precoding_settings)

    [ "Shi2011_WMMSE" => [ "user_rates" => user_rates1 ], 
      "Gomadam2008_MaxSINR" => [ "user_rates" => user_rates2 ],
      "Komulainen2013_WMMSE" => [ "user_rates" => user_rates3 ],
    ]
end

function plot_single_realization(results, simulation_params, precoding_settings)

    # Per-stream rate
    fig = plt.figure(figsize=(12,4*simulation_params["K"]))
    subplot_ind = 1

    for k = 1:simulation_params["K"]
        for n = 1:simulation_params["d"]
            ax = fig[:add_subplot](simulation_params["K"], simulation_params["d"], subplot_ind); subplot_ind += 1

            ax[:plot](squeeze(results["Shi2011_WMMSE"]["user_rates"][k,n,:], [1,2]), "b-")
            ax[:plot](squeeze(results["Gomadam2008_MaxSINR"]["user_rates"][k,n,:], [1,2]), "r-")
            ax[:plot](squeeze(results["Komulainen2013_WMMSE"]["user_rates"][k,n,:], [1,2]), "b-.")

            if n == 1
                ax[:set_ylabel]("User $k, rate")
            end
            if k == 1
                ax[:set_title]("Stream $n")
            end
            if k == simulation_params["K"]
                ax[:set_xlabel]("Time evolution")
            end
        end
    end
    fig[:savefig]("canonical_example_perstream.png", dpi=125)

    # Per-user sum rate
    fig = plt.figure(figsize=(12,4*simulation_params["K"]))
    subplot_ind = 1

    for k = 1:simulation_params["K"]
        ax = fig[:add_subplot](simulation_params["K"], 1, subplot_ind); subplot_ind += 1

        ax[:plot](squeeze(sum(results["Shi2011_WMMSE"]["user_rates"][k,:,:], 2), [1,2]), "b-")
        ax[:plot](squeeze(sum(results["Gomadam2008_MaxSINR"]["user_rates"][k,:,:], 2), [1,2]), "r-")
        ax[:plot](squeeze(sum(results["Komulainen2013_WMMSE"]["user_rates"][k,:,:], 2), [1,2]), "b-.")

        ax[:set_ylabel]("User $k, rate")

        if k == simulation_params["K"]
            ax[:set_xlabel]("Time evolution")
        end
    end
    fig[:savefig]("canonical_example_peruser.png", dpi=125)

    # System sum rate
    fig = plt.figure(figsize=(12,6))
    ax = fig[:add_subplot](1, 1, 1)

    ax[:plot](squeeze(sum(results["Shi2011_WMMSE"]["user_rates"], [1,2]), [1,2]), "b-", label="Shi2011_WMMSE")
    ax[:plot](squeeze(sum(results["Gomadam2008_MaxSINR"]["user_rates"], [1,2]), [1,2]), "r-", label="Gomadam2008_MaxSINR")
    ax[:plot](squeeze(sum(results["Komulainen2013_WMMSE"]["user_rates"], [1,2]), [1,2]), "b-.", label="Komulainen2013_WMMSE")

    ax[:set_xlabel]("SNR [dB]"); ax[:set_ylabel]("Sum rate [bits/s/Hz]")
    ax[:legend](loc="best")

    fig[:savefig]("canonical_example_system.png", dpi=125)
end

##########################################################################
# Settings
srand(23948189237)

const simulation_params = {
    "K" => 3, "N" => 3, "M" => 3,
    "P_dB" => 20, "sigma2_dB" => 0, "d" => 2
    }
const precoding_settings = {
    "stop_crit" => 50,
    "initial_precoders" => "dft"
    }

results = simulate_single_canonical_realization(simulation_params, precoding_settings)
plot_single_realization(results, simulation_params, precoding_settings)
