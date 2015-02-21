immutable Assignment
    cell_assignment::Vector{Int} # cell_assignment[k] is the BS serving MS k
    cell_assignment_inverse::Vector{IntSet} # cell_assignment_inverse[i] is the IntSet of MSs that BS i is serving

    cluster_assignment::Vector{IntSet} # cluster_assignment[k] is the IntSet of BSs coordinating to MS k
    cluster_assignment_inverse::Vector{IntSet} # cluster_assignment_inverse[i] is the IntSet of MSs that BS i is coordinating to
end
Assignment() = Assignment(Array(Int, 0), [IntSet()], [IntSet()], [IntSet()])

Assignment(cell_assignment::Vector{Int}, no_BSs::Int) =
    Assignment(cell_assignment, ones(length(cell_assignment), no_BSs))

Assignment(cell_assignment::Vector{Int}, cluster_assignment::Matrix) =
    Assignment(cell_assignment,
               [ IntSet(findin(cell_assignment, i)) for i = 1:size(cluster_assignment, 2) ],
               [ IntSet(findin(cluster_assignment[k,:], [1])) for k = 1:size(cluster_assignment, 1) ],
               [ IntSet(findin(cluster_assignment[:,i], [1])) for i = 1:size(cluster_assignment, 2) ]
    )

serving_BS_id(MS_id, assignment) =
    assignment.cell_assignment[MS_id]
served_MS_ids(BS_id, assignment) =
    assignment.cell_assignment_inverse[BS_id]
served_MS_ids_except_me(MS_id, BS_id, assignment) =
    setdiff(assignment.cell_assignment_inverse[BS_id], MS_id)

coordinated_BS_ids(MS_id, assignment) =
    assignment.cluster_assignment[MS_id]
coordinated_MS_ids(BS_id, assignment) =
    assignment.cluster_assignment_inverse[BS_id]

function require_equal_no_MSs_per_cell(assignment)
    BS_ind_max = maximum(cell_assignment.cell_assignment)

    Kc_test = length(served_MS_ids(1, cell_assignment))
    for BS_ind = 2:BS_ind_max
        (length(served_MS_ids(BS_ind, cell_assignment)) == Kc_test) || error("BSs must all serve equal number of MSs.")
    end
end
