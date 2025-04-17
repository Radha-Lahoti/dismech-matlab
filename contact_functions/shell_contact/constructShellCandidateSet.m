function [candidate_set, closest_dist] = constructShellCandidateSet(q, tri_combos, candidate_limit, face_nodes)
% Compute the min-distances of all possible edge combinations
minDs = min_distance_between_triangles(tri_combos, q);
closest_dist = min(minDs);

% Compute the indices of all edge combinations within the candidate_limit
col_indices = find(minDs < candidate_limit);
candidate_set = tri_combos(col_indices,:); % set C of edge_combos "in contact"

end

