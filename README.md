## Simulation in Wireless Communications: CoordinatedPrecoding.jl
CoordinatedPrecoding.jl is [Julia][julia] package for the simulation of
wireless networks, in particular multicell multiple-input multiple-output (MIMO)
networks operating using [coordinated precoding][PracticalTDD].

The main purpose of the package is to provide an infrastructure for researchers
to write their own simulation scripts, for example to try out new ideas for
coordinated precoding. In addition to this, the package contains some basic
network/channel models (inspired by 3GPP and ITU-T) and some reference
implementations of popular precoding methods.

### Overview

#### Networks
A network to be studied is given by `network::Network`. In order to generate
the network, several models are available:

* `InterferingBroadcastChannel`: network without large-scale fading which
  can be used to simulate interfering broadcast channels
  (including interference channels).
* `TriangularMacroNetwork`: a triangular macro network with three sites.
  Large-scale parameters are inspired by 3GPP Case 1.
* `TriangularHetNet`: similar to the network above, but including a variable
  number of pico base stations per cell.
* `IndoorsNetwork`: an office corridor inspired by ITU-T InH
* `RandomLargeScaleNetwork`: a macro network with random base station placement
  and 3GPP inspired pathloss.

For each network, a convenience function is available for quick generation, e.g.:
`setup_interfering_broadcast_channel(num_BSs, num_MSs_per_cell, num_MS_antennas, num_BS_antennas; kwargs)`

#### Simulation functions
The simulation is a simple Monte Carlo loop. Given a network, the
`simulate(network, simulation_params)` function determines large-scale
parameters by dropping the users by calling `draw_user_drop!(network)`. The
small-scale fading is then generated by calling `draw_channel(network)`. By
writing customized methods for these functions, new network models can be
implemented.

The simulation function sweeps one main independent variable, but support
also exists for sweeping auxiliary independent variables.

### Reference implementations of precoding schemes
We provide implementations for the following precoding methods:

* `Shi2011_WMMSE`: the famous WMMSE algorithm of [Shi et al.][Shi2011]
* `Razaviyayn2013_MinMaxWMMSE`: min-max version of the WMMSE algorithm by [Razaviyayn et al.][Razaviyayn2013]
* `Komulainen2013_WMMSE`: diagonalized WMMSE algorithm of [Komulainen et al.][Komulainen2013]
* `Gomadam2008_MinWLI`: leakage minimization of [Gomadam et al.][Gomadam2008]
* `Gomadam2008_MaxSINR`: heuristic per-stream SINR maximization of [Gomadam et al.][Gomadam2008]
* `Eigenprecoding`: single-user optimal beamforming with/without TDMA

### Usage
See the [examples](examples) folder for ideas on how to write the simulation
scripts. For implementing network models and precoding methods, see the
corresponding folders in the `src` directory.

The package uses PyPlot.jl for visualization, Lumberjack for logging, and
ProgressMeter for showing progress during lengthy simulation runs. The Min-Max
WMMSE algorithm is implemented using Gurobi.

### Improvements
The package is under continuous development. Several aspects of the code could
be improved, e.g. test coverage and harmonization of concepts. This is not
likely to happen in the short term however.

### License
This source code is licensed under the X license.

[julia]: http://www.julialang.org
[PracticalTDD]: http://kth.diva-portal.org/smash/get/diva2:811008/FULLTEXT01.pdf
[Shi2011]: http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=5756489
[Razaviyayn2013]: http://www.sciencedirect.com/science/article/pii/S0165168413000716
[Komulainen2013]: http://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=6463462
[Gomadam2008]: http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=5773023
