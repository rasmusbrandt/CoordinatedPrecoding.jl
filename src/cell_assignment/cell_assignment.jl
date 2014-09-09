immutable CellAssignment
    assignment::Vector{Int}
end

get_cell_assignment(network::Network) = CellAssignment(Int[ findin(network.BSs, [network.MSs[k].served_by_BS])[1] for k = 1:get_no_MSs(network) ])

serving_BS_id(MS_id::Int, cell_assignment::CellAssignment) = cell_assignment.assignment[MS_id]
served_MS_ids(BS_id::Int, cell_assignment::CellAssignment) = IntSet(findin(cell_assignment.assignment, BS_id))

function served_MS_ids_except_me(MS_id::Int, BS_id::Int, cell_assignment::CellAssignment)
    all_served_MS_ids = served_MS_ids(BS_id, cell_assignment)
    delete!(all_served_MS_ids, MS_id)
end

function assign_cells_by_pathloss!(network::PhysicalNetwork)
    # Sets cell assignment (in the MS objects) based on pathloss of physical network
    error("Implement me")
end

function assign_cells_by_instantaneous_channels!{Channel_t <: Channel}(network::PhysicalNetwork, channel::Channel_t)
    # Sets cell assignment (in the MS objects) based on instantaneous channel strength
    error("Implement me")
end
