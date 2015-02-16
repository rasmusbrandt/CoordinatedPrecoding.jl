immutable ClusterAssignment
    assignment::Vector{IntSet} # assignment[k] is the IntSet of BSs coordinating to MS k
    assignment_inverse::Vector{IntSet} # assignment_inverse[i] is the IntSet of MSs that BS i is coordinating to
end
ClusterAssignment() = ClusterAssignment([IntSet()], [IntSet()])
ClusterAssignment(assignment::Matrix) =
    ClusterAssignment([ IntSet(findin(assignment[k,:], [1])) for k = 1:size(assignment, 1) ],
                      [ IntSet(findin(assignment[:,l], [1])) for l = 1:size(assignment, 2) ])

coordinated_BS_ids(MS_id::Int, cluster_assignment::ClusterAssignment) =
    cluster_assignment.assignment[MS_id]
coordinated_MS_ids(BS_id::Int, cluster_assignment::ClusterAssignment) =
    cluster_assignment.assignment_inverse[BS_id]
