function [dist, is_inside] = pointToTriangleDist(p, tri)
    a = tri(:,1);
    b = tri(:,2);
    c = tri(:,3);
    n = cross(b - a, c - a);
    n = n / norm(n);
    dist = dot(p - a, n);
    proj = p - dist * n;
    is_inside = pointInTriangle(proj, tri);
    dist = abs(dist);
end

function inside = pointInTriangle(p, tri)
    a = tri(:,1); b = tri(:,2); c = tri(:,3);
    v0 = c - a;
    v1 = b - a;
    v2 = p - a;
    d00 = dot(v0, v0);
    d01 = dot(v0, v1);
    d02 = dot(v0, v2);
    d11 = dot(v1, v1);
    d12 = dot(v1, v2);
    invDenom = 1 / (d00 * d11 - d01 * d01);
    u = (d11 * d02 - d01 * d12) * invDenom;
    v = (d00 * d12 - d01 * d02) * invDenom;
    inside = (u >= 0) && (v >= 0) && (u + v <= 1);
end