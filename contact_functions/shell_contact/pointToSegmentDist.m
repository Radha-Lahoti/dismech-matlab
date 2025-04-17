function [dist, closest_point] = pointToSegmentDist(p, a, b)
    ab = b - a;
    t = dot(p - a, ab) / dot(ab, ab);
    t = max(0, min(1, t));
    closest_point = a + t * ab;
    dist = norm(p - closest_point);
end
