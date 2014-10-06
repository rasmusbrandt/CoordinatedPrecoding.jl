using CoordinatedPrecoding
using PyCall; @pyimport matplotlib.pyplot as plt

##########################################################################
# Simulation and visualization functions
function simulate_single(simulation_params, precoding_settings)
    network = setup_interference_channel(simulation_params["K"],
        simulation_params["N"], simulation_params["M"],
        transmit_power=10^(simulation_params["P_dB"]/10),
        receiver_noise_power=10^(simulation_params["sigma2_dB"]/10),
        no_streams=simulation_params["d"])

    cell_assignment = assign_cells_by_id(network)

    channel = draw_channel(network)

    user_rates1 = Shi2011_WMMSE(channel, network, cell_assignment, precoding_settings)
    user_rates2 = Gomadam2008_MaxSINR(channel, network, cell_assignment, precoding_settings)
    user_rates3 = Komulainen2013_WMMSE(channel, network, cell_assignment, precoding_settings)
    user_rates_uncoord, user_rates_intercell_tdma, user_rates_intracell_tdma = Eigenprecoding(channel, network, cell_assignment, precoding_settings)

    [ "Shi2011_WMMSE" => [ "user_rates" => user_rates1 ], 
      "Gomadam2008_MaxSINR" => [ "user_rates" => user_rates2 ],
      "Komulainen2013_WMMSE" => [ "user_rates" => user_rates3 ],
      "Eigenprecoding" => [ "user_rates_uncoord" => user_rates_uncoord,
                            "user_rates_intercell_tdma" => user_rates_intercell_tdma,
                            "user_rates_intracell_tdma" => user_rates_intracell_tdma ],
    ]
end

function plot_single(results, simulation_params, precoding_settings)
    iters = 1:precoding_settings["stop_crit"]

    # Per-stream rate
    fig = plt.figure(figsize=(12,4*simulation_params["K"]))
    subplot_ind = 1

    for k = 1:simulation_params["K"]
        for n = 1:simulation_params["d"]
            ax = fig[:add_subplot](simulation_params["K"], simulation_params["d"], subplot_ind); subplot_ind += 1

            ax[:plot](iters, squeeze(results["Shi2011_WMMSE"]["user_rates"][k,n,:], [1,2]), "b-")
            ax[:plot](iters, squeeze(results["Gomadam2008_MaxSINR"]["user_rates"][k,n,:], [1,2]), "r-")
            ax[:plot](iters, squeeze(results["Komulainen2013_WMMSE"]["user_rates"][k,n,:], [1,2]), "b-.")
            ax[:plot](iters, results["Eigenprecoding"]["user_rates_intercell_tdma"][k,n]*ones(iters), "m-")
            ax[:plot](iters, results["Eigenprecoding"]["user_rates_intracell_tdma"][k,n]*ones(iters), "m.-")
            ax[:plot](iters, results["Eigenprecoding"]["user_rates_uncoord"][k,n]*ones(iters), "k-")

            if n == 1
                ax[:set_ylabel]("User $k")
            end
            if k == 1
                ax[:set_title]("Stream $n")
            end
            if k == simulation_params["K"]
                ax[:set_xlabel]("Time evolution")
            end
        end
    end
    fig[:savefig]("interference_channel_single_perstream.png", dpi=125)

    # Per-user sum rate
    fig = plt.figure(figsize=(12,4*simulation_params["K"]))
    subplot_ind = 1

    for k = 1:simulation_params["K"]
        ax = fig[:add_subplot](simulation_params["K"], 1, subplot_ind); subplot_ind += 1

        ax[:plot](iters, squeeze(sum(results["Shi2011_WMMSE"]["user_rates"][k,:,:], 2), [1,2]), "b-")
        ax[:plot](iters, squeeze(sum(results["Gomadam2008_MaxSINR"]["user_rates"][k,:,:], 2), [1,2]), "r-")
        ax[:plot](iters, squeeze(sum(results["Komulainen2013_WMMSE"]["user_rates"][k,:,:], 2), [1,2]), "b-.")
        ax[:plot](iters, sum(results["Eigenprecoding"]["user_rates_intercell_tdma"][k,:])*ones(iters), "m-")
        ax[:plot](iters, sum(results["Eigenprecoding"]["user_rates_intracell_tdma"][k,:])*ones(iters), "m.-")
        ax[:plot](iters, sum(results["Eigenprecoding"]["user_rates_uncoord"][k,:])*ones(iters), "k-")

        ax[:set_ylabel]("User $k")

        if k == simulation_params["K"]
            ax[:set_xlabel]("Time evolution")
        end
    end
    fig[:savefig]("interference_channel_single_peruser.png", dpi=125)

    # System sum rate
    fig = plt.figure(figsize=(12,6))
    ax = fig[:add_subplot](1, 1, 1)

    ax[:plot](iters, squeeze(sum(results["Shi2011_WMMSE"]["user_rates"], [1,2]), [1,2]), "b-", label="Shi2011_WMMSE")
    ax[:plot](iters, squeeze(sum(results["Gomadam2008_MaxSINR"]["user_rates"], [1,2]), [1,2]), "r-", label="Gomadam2008_MaxSINR")
    ax[:plot](iters, squeeze(sum(results["Komulainen2013_WMMSE"]["user_rates"], [1,2]), [1,2]), "b-.", label="Komulainen2013_WMMSE")
    ax[:plot](iters, sum(results["Eigenprecoding"]["user_rates_intercell_tdma"])*ones(iters), "m-", label="Intercell TDMA")
    ax[:plot](iters, sum(results["Eigenprecoding"]["user_rates_intracell_tdma"])*ones(iters), "m.-", label="Intracell TDMA")
    ax[:plot](iters, sum(results["Eigenprecoding"]["user_rates_uncoord"])*ones(iters), "k-", label="Uncoordinated")

    ax[:set_xlabel]("Time evolution"); ax[:set_ylabel]("Sum rate [bits/s/Hz]")
    ax[:legend](loc="best")

    fig[:savefig]("interference_channel_single_system.png", dpi=125)
end

##########################################################################
# Settings
srand(int(time()))

const simulation_params = {
    "K" => 3, "N" => 3, "M" => 3,
    "P_dB" => 20, "sigma2_dB" => 0, "d" => 2
}
const precoding_settings = {
    "stop_crit" => 50,
    "initial_precoders" => "dft"
}

results = simulate_single(simulation_params, precoding_settings)
plot_single(results, simulation_params, precoding_settings)
