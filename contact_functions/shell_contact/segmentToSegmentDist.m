function [dist, p1, p2, ratio1, ratio2] = segmentToSegmentDist(a1, b1, a2, b2)
    u = b1 - a1;
    v = b2 - a2;
    w0 = a1 - a2;

    a = dot(u, u);
    b = dot(u, v);
    c = dot(v, v);
    d = dot(u, w0);
    e = dot(v, w0);

    D = a * c - b * b;
    s = 0;
    t = 0;

    if D > 1e-8
        s = (b * e - c * d) / D;
        t = (a * e - b * d) / D;
        s = max(0, min(1, s));
        t = max(0, min(1, t));
    end

    p1 = a1 + s * u;
    p2 = a2 + t * v;
    dist = norm(p1 - p2);

    ratio1 = [1-s, s]; %[a1, b1]
    ratio2 = [1-t, t]; %[a2, b2]
end