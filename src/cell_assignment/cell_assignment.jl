immutable CellAssignment
    assignment::Vector{Int} # length is number of MSs
    assignment_inverse::Vector{IntSet} # length is number of BSs
end
CellAssignment(assignment::Vector{Int}, no_BSs::Int) =
    CellAssignment(assignment, [ IntSet(findin(assignment, BS_id)) for BS_id = 1:no_BSs ])

serving_BS_id(MS_id::Int, cell_assignment::CellAssignment) =
    cell_assignment.assignment[MS_id]
served_MS_ids(BS_id::Int, cell_assignment::CellAssignment) =
    cell_assignment.assignment_inverse[BS_id]
served_MS_ids_except_me(MS_id::Int, BS_id::Int, cell_assignment::CellAssignment) =
    setdiff(cell_assignment.assignment_inverse[BS_id], MS_id)

function assign_cells_by_pathloss!(network::PhysicalNetwork)
    # Cell assignment based on pathloss of physical network
    error("Implement me")
end

function assign_cells_by_instantaneous_channels!{Channel_t <: Channel}(network::PhysicalNetwork, channel::Channel_t)
    # Cell assignment based on instantaneous channel strength
    error("Implement me")
end

function require_equal_no_MSs_per_cell(cell_assignment::CellAssignment)
    BS_ind_max = maximum(cell_assignment.assignment)

    Kc_test = length(served_MS_ids(1, cell_assignment))
    for BS_ind = 2:BS_ind_max
        (length(served_MS_ids(BS_ind, cell_assignment)) == Kc_test) || error("BSs must all serve equal number of MSs.")
    end
end
