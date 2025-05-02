function Econ = computeShellContactEnergy(data, constraint_type, h, delta, K1)
% Extract position data
x11 = data(1:3);
x12 = data(4:6);
x13 = data(7:9);
x21 = data(10:12);
x22 = data(13:15);
x23 = data(16:18);

if (constraint_type=="PointToPoint")
    Delta = norm(x21-x11);

elseif (constraint_type=="PointToEdge")
    % point to edge
    e2 = x22 - x21; % the edge for which the minimum distance vector does not lie on
    % an end
    e_other = x21-x11; % x11 is the node of edge 1 which the minimum distance vector does lie on
    Delta = norm(cross(e2,e_other))/norm(e2);

elseif (constraint_type=="EdgeToEdge")
    % edge to edge
    e1 = x12 - x11;
    e2 = x22 - x21;
    u = cross(e1,e2);
    u_hat = u/norm(u);
    Delta = abs(dot((x21-x11),u_hat));

elseif (constraint_type=="PointToFace")
    % point to triangle 
    % distance from vertex on triangle 1 to some point inside triangle 2
    n2 = cross((x22-x21),(x23-x21));
    n2_hat = n2/norm(n2); % normal of triangle 2
    Delta = abs(dot((x21-x11),n2_hat));

elseif (constraint_type=="EdgeToFace")
    Econ = []; % random
    return;
end

% Define Energy 
if Delta <= 2*h - delta
    Econ = (2*h-Delta)^2;
elseif Delta > 2*h - delta && Delta < 2*h + delta
    Econ = ((1/K1)*log(1+exp(K1*(2*h-Delta))))^2;
end
