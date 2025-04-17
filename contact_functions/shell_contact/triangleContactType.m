function [dist, contact_type, tri_pair_updated, tri_pair_nodes_updated] = triangleContactType(q, tri_pair_nodes, scale)

% Returns:
%   dist: minimum distance between the two triangles
%   contact_type: one of {'PointToPoint', 'PointToEdge', 'EdgeToEdge', 'PointToFace', 'EdgeToFace'}
%   contact_info: structure with relevant points or segments
%   tri_pair_updated: rearranged triangle vertices

    % All vertices
    face_nodes1 = tri_pair_nodes(1:3);
    face_nodes2 = tri_pair_nodes(4:6);
    V1 = scale.*[q(mapNodetoDOF(face_nodes1(1))) , q(mapNodetoDOF(face_nodes1(2))) , q(mapNodetoDOF(face_nodes1(3)))];
    V2 = scale.*[q(mapNodetoDOF(face_nodes2(1))) , q(mapNodetoDOF(face_nodes2(2))) , q(mapNodetoDOF(face_nodes2(3)))];
    % V1 and V2 are 3x3 matrices: each column is a vertex of the triangle in 3D

    min_dist = inf;
    contact_type = "";
    contact_info = struct();

    tri_pair_nodes_updated = [face_nodes1, face_nodes2];
    tri_pair_updated = [reshape(V1,[9,1]); reshape(V2,[9,1])];

    updated_tri1 = V1;
    updated_tri2 = V2;
    updated_face_nodes1 = face_nodes1;
    updated_face_nodes2 = face_nodes2;
    idx1_pp = -1;
    idx2_pp = -1;

    % --- Point-to-Point distances ---
    for i = 1:3
        for j = 1:3
            d = norm(V1(:,i) - V2(:,j));
            if d < min_dist
                min_dist = d;
                contact_type = "PointToPoint";
                contact_info = struct('p1', V1(:,i), 'p2', V2(:,j));
                idx1_pp = i;
                idx2_pp = j;
            end
        end
    end

    % --- Point-to-Edge distances ---
    for i = 1:3
        for j = 1:3
            a = V2(:,j);
            b = V2(:,mod(j,3)+1);
            [d, cp] = pointToSegmentDist(V1(:,i), a, b);
            if d < min_dist
                min_dist = d;
                contact_type = "PointToEdge";
                contact_info = struct('point', V1(:,i), 'edge', [a, b], 'closest_point', cp);
            end
        end
    end

    for i = 1:3
        for j = 1:3
            a = V1(:,j);
            b = V1(:,mod(j,3)+1);
            [d, cp] = pointToSegmentDist(V2(:,i), a, b);
            if d < min_dist
                min_dist = d;
                contact_type = "PointToEdge";
                contact_info = struct('point', V2(:,i), 'edge', [a, b], 'closest_point', cp);
            end
        end
    end

    % --- Edge-to-Edge distances ---
    for i = 1:3
        a1 = V1(:,i);
        b1 = V1(:,mod(i,3)+1);
        for j = 1:3
            a2 = V2(:,j);
            b2 = V2(:,mod(j,3)+1);
            [d, p1, p2] = segmentToSegmentDist(a1, b1, a2, b2);
            if d < min_dist
                min_dist = d;
                contact_type = "EdgeToEdge";
                contact_info = struct('edge1', [a1, b1], 'edge2', [a2, b2], 'p1', p1, 'p2', p2);
            end
        end
    end

    % --- Point-to-Face (Point inside Triangle) ---
    for i = 1:3
        [d, is_inside] = pointToTriangleDist(V1(:,i), V2);
        if is_inside && d < min_dist
            min_dist = d;
            contact_type = "PointToFace";
            contact_info = struct('point', V1(:,i), 'triangle', V2);
        end

        [d, is_inside] = pointToTriangleDist(V2(:,i), V1);
        if is_inside && d < min_dist
            min_dist = d;
            contact_type = "PointToFace";
            contact_info = struct('point', V2(:,i), 'triangle', V1);
        end
    end

    % --- Edge-to-Face (Intersection case) ---
    for i = 1:3
        a = V1(:,i);
        b = V1(:,mod(i,3)+1);
        [intersects, pt] = segmentTriangleIntersection(a, b, V2);
        if intersects
            if norm(pt - a) > 1e-8 && norm(pt - b) > 1e-8
                min_dist = 0;
                contact_type = "EdgeToFace";
                contact_info = struct('edge', [a, b], 'triangle', V2, 'intersection', pt);
                dist = min_dist;
                return
            end
        end

        a = V2(:,i);
        b = V2(:,mod(i,3)+1);
        [intersects, pt] = segmentTriangleIntersection(a, b, V1);
        if intersects
            if norm(pt - a) > 1e-8 && norm(pt - b) > 1e-8
                min_dist = 0;
                contact_type = "EdgeToFace";
                contact_info = struct('edge', [a, b], 'triangle', V1, 'intersection', pt);
                dist = min_dist;
                return
            end
        end
    end

    dist = min_dist;

    % Rearranging triangles and face nodes for contact types
    if contact_type == "PointToPoint"
        updated_tri1 = circshift(V1, [0, 1 - idx1_pp]);
        updated_face_nodes1 = circshift(face_nodes1, [0, 1 - idx1_pp]);
        updated_tri2 = circshift(V2, [0, 1 - idx2_pp]);
        updated_face_nodes2 = circshift(face_nodes2, [0, 1 - idx2_pp]);

    elseif contact_type == "PointToEdge"
        point = contact_info.point;
        edge = contact_info.edge;

        if ismember(point', V1', 'rows')
            idx_point = find(ismember(V1', point', 'rows'));
            updated_tri1 = circshift(V1, [0, 1 - idx_point]);
            updated_face_nodes1 = circshift(face_nodes1, [0, 1 - idx_point]);
            idx_edge_start = find(ismember(V2', edge(:,1)', 'rows'));
            updated_tri2 = circshift(V2, [0, 1 - idx_edge_start]);
            updated_face_nodes2 = circshift(face_nodes2, [0, 1 - idx_edge_start]);
        else
            idx_point = find(ismember(V2', point', 'rows'));
            updated_tri2 = circshift(V2, [0, 1 - idx_point]);
            updated_face_nodes2 = circshift(face_nodes2, [0, 1 - idx_point]);
            idx_edge_start = find(ismember(V1', edge(:,1)', 'rows'));
            updated_tri1 = circshift(V1, [0, 1 - idx_edge_start]);
            updated_face_nodes1 = circshift(face_nodes1, [0, 1 - idx_edge_start]);
        end

    elseif contact_type == "EdgeToEdge"
        edge1 = contact_info.edge1;
        edge2 = contact_info.edge2;
        idx1 = find(ismember(V1', edge1(:,1)', 'rows'));
        updated_tri1 = circshift(V1, [0, 1 - idx1]);
        updated_face_nodes1 = circshift(face_nodes1, [0, 1 - idx1]);
        idx2 = find(ismember(V2', edge2(:,1)', 'rows'));
        updated_tri2 = circshift(V2, [0, 1 - idx2]);
        updated_face_nodes2 = circshift(face_nodes2, [0, 1 - idx2]);

    elseif contact_type == "PointToFace"
        point = contact_info.point;
        if ismember(point', V1', 'rows')
            idx_point = find(ismember(V1', point', 'rows'));
            updated_tri1 = circshift(V1, [0, 1 - idx_point]);
            updated_face_nodes1 = circshift(face_nodes1, [0, 1 - idx_point]);
        else
            idx_point = find(ismember(V2', point', 'rows'));
            updated_tri2 = circshift(V2, [0, 1 - idx_point]);
            updated_face_nodes2 = circshift(face_nodes2, [0, 1 - idx_point]);
        end

    elseif contact_type == "EdgeToFace"
        edge = contact_info.edge;
        if ismember(edge(:,1)', V1', 'rows')
            idx_edge_start = find(ismember(V1', edge(:,1)', 'rows'));
            updated_tri1 = circshift(V1, [0, 1 - idx_edge_start]);
            updated_face_nodes1 = circshift(face_nodes1, [0, 1 - idx_edge_start]);
        else
            idx_edge_start = find(ismember(V2', edge(:,1)', 'rows'));
            updated_tri2 = circshift(V2, [0, 1 - idx_edge_start]);
            updated_face_nodes2 = circshift(face_nodes2, [0, 1 - idx_edge_start]);
        end
    end

    tri_pair_updated = [reshape(updated_tri1,[9,1]); reshape(updated_tri2,[9,1])];
    tri_pair_nodes_updated = [updated_face_nodes1, updated_face_nodes2];
end
