function [in_contact_tri_combos, in_contact_constraints, closest_dist, num_coll] = detectShellCollisions_trial(q, tri_combos,  constraints, delta, contact_len)
% Compute the min-distances of all possible edge combinations
minDs = min_distance_between_triangles(tri_combos, q);
closest_dist = min(minDs); 

% Compute the indices of all edge combinations within the collision limit
col_indices = find(minDs < (delta + contact_len));
in_contact_tri_combos = tri_combos(col_indices,:); % set C of edge_combos "in contact"
in_contact_constraints = constraints(col_indices);

num_coll = size(in_contact_tri_combos); % number of collisions

end

