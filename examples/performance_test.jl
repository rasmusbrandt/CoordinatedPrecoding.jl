using CoordinatedPrecoding

##########################################################################
# Performance functions
function run_performance(simulation_params, precoding_settings, performance_params)
    network = setup_interference_channel(simulation_params["K"],
        simulation_params["N"], simulation_params["M"],
        transmit_power=10^(simulation_params["P_dB"]/10),
        receiver_noise_power=10^(simulation_params["sigma2_dB"]/10),
        no_streams=simulation_params["d"])

    channel = draw_channel(network)

    Shi2011_WMMSE(channel, network, precoding_settings)
    @time for i = 1:performance_params["Nperf"]; Shi2011_WMMSE(channel, network, precoding_settings); end
    Gomadam2008_MaxSINR(channel, network, precoding_settings)
    @time for i = 1:performance_params["Nperf"]; Gomadam2008_MaxSINR(channel, network, precoding_settings); end
    Komulainen2013_WMMSE(channel, network, precoding_settings)
    @time for i = 1:performance_params["Nperf"]; Komulainen2013_WMMSE(channel, network, precoding_settings); end
end


##########################################################################
# Settings
const simulation_params = {
    "K" => 3, "N" => 2, "M" => 2,
    "P_dB" => 30, "sigma2_dB" => 0, "d" => 1
}
const precoding_settings = { "stop_crit" => 50 }
const performance_params = { "Nperf" => 1000 }

run_performance(simulation_params, precoding_settings, performance_params)
