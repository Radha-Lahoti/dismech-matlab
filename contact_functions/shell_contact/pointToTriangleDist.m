function [dist, is_inside, baryCoords] = pointToTriangleDist(p, tri)
    a = tri(:,1);
    b = tri(:,2);
    c = tri(:,3);
    
    % Triangle normal
    n = cross(b - a, c - a);
    n = n / norm(n);
    
    % Project point onto triangle plane
    dist = dot(p - a, n);
    proj = p - dist * n;
    dist = abs(dist);
    
    % Compute barycentric coordinates
    v0 = b - a;
    v1 = c - a;
    v2 = proj - a;
    
    d00 = dot(v0, v0);
    d01 = dot(v0, v1);
    d11 = dot(v1, v1);
    d20 = dot(v2, v0);
    d21 = dot(v2, v1);
    
    denom = d00 * d11 - d01 * d01;
    v = (d11 * d20 - d01 * d21) / denom;
    w = (d00 * d21 - d01 * d20) / denom;
    u = 1 - v - w;
    
    baryCoords = [u; v; w];  % Corresponds to weights for a, b, c
    is_inside = all(baryCoords >= 0) && all(baryCoords <= 1);
end
