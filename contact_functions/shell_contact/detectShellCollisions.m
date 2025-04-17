function [in_contact_tri_combos, closest_dist, num_coll] = detectShellCollisions(q, tri_combos,  delta, contact_len, face_nodes)
% Compute the min-distances of all possible edge combinations
minDs = min_distance_between_triangles(tri_combos, q);
closest_dist = min(minDs); 

% Compute the indices of all edge combinations within the collision limit
col_indices = find(minDs < (delta + contact_len));
in_contact_tri_combos = tri_combos(col_indices,:); % set C of edge_combos "in contact"

num_coll = size(in_contact_tri_combos); % number of collisions

end

