function [intersects, pt] = segmentTriangleIntersection(p0, p1, tri)
    a = tri(:,1); b = tri(:,2); c = tri(:,3);
    dir = p1 - p0;
    edge1 = b - a;
    edge2 = c - a;
    h = cross(dir, edge2);
    det = dot(edge1, h);
    if abs(det) < 1e-8
        intersects = false;
        pt = [];
        return;
    end
    inv_det = 1 / det;
    s = p0 - a;
    u = dot(s, h) * inv_det;
    q = cross(s, edge1);
    v = dot(dir, q) * inv_det;
    t = dot(edge2, q) * inv_det;
    intersects = (u >= 0) && (v >= 0) && (u + v <= 1) && (t >= 0) && (t <= 1);
    pt = p0 + t * dir;
end