function minDs = min_distance_between_triangles(tri_combos, q)

assert(size(tri_combos, 2) == 6);
num_inputs = size(tri_combos, 1);
tri_combo_inputs = zeros(size(tri_combos,1),18);

for i = 1:size(tri_combos,1)
    for j = 1:size(tri_combos,2)
        tri_combo_inputs(i,3*j-2:3*j) = q(mapNodetoDOF(tri_combos(i,j)));
    end
end

min_dist = inf;
minDs = zeros(num_inputs,1);

for c = 1:num_inputs

    % All vertices
    V1 = [ q(mapNodetoDOF(tri_combos(c, 1))) , q(mapNodetoDOF(tri_combos(c, 2))) , q(mapNodetoDOF(tri_combos(c, 3))) ];
    V2 = [ q(mapNodetoDOF(tri_combos(c, 4))) , q(mapNodetoDOF(tri_combos(c, 5))) , q(mapNodetoDOF(tri_combos(c, 6))) ];


    % --- Point-to-Point distances ---
    for i = 1:3
        for j = 1:3
            d = norm(V1(:,i) - V2(:,j));
            if d < min_dist
                min_dist = d;
            end
        end
    end

    % --- Point-to-Edge distances ---
    for i = 1:3
        for j = 1:3
            a = V2(:,j);
            b = V2(:,mod(j,3)+1);
            [d] = pointToSegmentDist(V1(:,i), a, b);
            if d < min_dist
                min_dist = d;
            end
        end
    end

    for i = 1:3
        for j = 1:3
            a = V1(:,j);
            b = V1(:,mod(j,3)+1);
            [d] = pointToSegmentDist(V2(:,i), a, b);
            if d < min_dist
                min_dist = d;
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
            [d] = segmentToSegmentDist(a1, b1, a2, b2);
            if d < min_dist
                min_dist = d;
            end
        end
    end

    % --- Point-to-Face (Point projection inside Triangle) ---
    for i = 1:3
        [d, is_inside] = pointToTriangleDist(V1(:,i), V2);
        if is_inside && d < min_dist
            min_dist = d;
        end

        [d, is_inside] = pointToTriangleDist(V2(:,i), V1);
        if is_inside && d < min_dist
            min_dist = d;
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
                minDs(c)= min_dist;
                break
            end
        end

        a = V2(:,i);
        b = V2(:,mod(i,3)+1);
        [intersects, pt] = segmentTriangleIntersection(a, b, V1);
        if intersects
            if norm(pt - a) > 1e-8 && norm(pt - b) > 1e-8
                min_dist = 0;
                minDs(c)= min_dist;
                break
            end
        end
    end

    minDs(c)= min_dist;
end