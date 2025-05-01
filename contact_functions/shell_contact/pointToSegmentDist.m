function [dist, closest_point, ratios] = pointToSegmentDist(p, a, b)
    ab = b - a;
    t = dot(p - a, ab) / dot(ab, ab);
    t = max(0, min(1, t));
    closest_point = a + t * ab;
    dist = norm(p - closest_point);
    ratios = [norm(b-closest_point)/norm(ab), norm(a-closest_point)/norm(ab)]; 
end
