function [grad_E1_p2p, grad_E1_p2e, grad_E1_e2e, grad_E1_p2t, grad_E2_p2p, grad_E2_p2e, grad_E2_e2e, grad_E2_p2t ]...
    = generateContactPiecewiseFunctions()

    % Declare symbolic variables
    syms x11 [3 1] real
    syms x12 [3 1] real
    syms x13 [3 1] real
    syms x21 [3 1] real
    syms x22 [3 1] real
    syms x23 [3 1] real
    syms Delta real
    syms h K1 real

    % point to point
    Delta_p2p = norm(x21-x11); % x11 and x21 are points of triangle 1 and 2 in contact
    % we need to make sure that x11 and x21 are correctly put during simulation

    % point to edge
    e2 = x22 - x21; % the edge for which the minimum distance vector does not lie on
    % an end
    e_other = x21-x11; % x11 is the node of edge 1 which the minimum distance vector does lie on
    Delta_p2e = norm(cross(e2,e_other))/norm(e2);

    % edge to edge
    e1 = x12 - x11;
    u = cross(e1,e2);
    u_hat = u/norm(u);
    Delta_e2e = abs(dot((x21-x11),u_hat));

    % point to triangle 
    % distance from vertex on triangle 1 to some point inside triangle 2
    n2 = cross((x22-x21),(x23-x21));
    n2_hat = n2/norm(n2); % normal of triangle 2
    Delta_p2t = abs(dot((x21-x11),n2_hat));

    % find gradient of Delta wrt [x11; x12; x13; x21; x22; x23]
    x = [x11; x12; x13; x21; x22; x23];
    % grad_Delta_p2p = (jacobian(Delta_p2p, x));
    grad_Delta_p2p = [x11./Delta_p2p ; zeros(3,1); zeros(3,1); x21./Delta_p2p; zeros(3,1); zeros(3,1)];

    grad_Delta_p2e = (jacobian(Delta_p2e, x));
    grad_Delta_e2e = (jacobian(Delta_e2e, x));
    grad_Delta_p2t = (jacobian(Delta_p2t, x));

    % find hessian of Delta wrt [x11; x12; x13; x21; x22; x23]
    hess_Delta_p2p = jacobian(grad_Delta_p2p, x);
    hess_Delta_p2e = jacobian(grad_Delta_p2e, x);
    hess_Delta_e2e = jacobian(grad_Delta_e2e, x);
    hess_Delta_p2t = jacobian(grad_Delta_p2t, x);

    % Define Energy 
    % if Δ <= 2h - δ
    E1 = (2*h-Delta)^2;
    % if (2h - δ) < Δ < (2h + δ)
    E2 = ((1/K1)*log(1+exp(K1*(2*h-Delta))))^2;

    % find gradient of E
    del_E1_del_Delta = diff(E1,Delta);
    del_E2_del_Delta = diff(E2,Delta);

    % find second derivative of E wrt Delta
    del2_E1_del_Delta2 = diff(E1,Delta,2);
    del2_E2_del_Delta2 = diff(E2,Delta,2);

    grad_E1_p2p = del_E1_del_Delta .* grad_Delta_p2p;
    grad_E1_p2e = del_E1_del_Delta .* grad_Delta_p2e;
    grad_E1_e2e = del_E1_del_Delta .* grad_Delta_e2e;
    grad_E1_p2t = del_E1_del_Delta .* grad_Delta_p2t;

    grad_E2_p2p = del_E2_del_Delta .* grad_Delta_p2p;
    grad_E2_p2e = del_E2_del_Delta .* grad_Delta_p2e;
    grad_E2_e2e = del_E2_del_Delta .* grad_Delta_e2e;
    grad_E2_p2t = del_E2_del_Delta .* grad_Delta_p2t;

    % hessian of E
    hess_E1_p2p = del_E1_del_Delta .* hess_Delta_p2p + ...
        (transpose(grad_Delta_p2p)*grad_Delta_p2p).*del2_E1_del_Delta2;
    hess_E1_p2e = del_E1_del_Delta .* hess_Delta_p2e + ...
        (transpose(grad_Delta_p2e)*grad_Delta_p2e).*del2_E1_del_Delta2;
    hess_E1_e2e = del_E1_del_Delta .* hess_Delta_e2e + ...
        (transpose(grad_Delta_e2e)*grad_Delta_e2e).*del2_E1_del_Delta2;
    hess_E1_p2t = del_E1_del_Delta .* hess_Delta_p2t + ...
        (transpose(grad_Delta_p2t)*grad_Delta_p2t).*del2_E1_del_Delta2;

    hess_E2_p2p = del_E2_del_Delta .* hess_Delta_p2p + ...
        (transpose(grad_Delta_p2p)*grad_Delta_p2p).*del2_E2_del_Delta2;
    hess_E2_p2e = del_E2_del_Delta .* hess_Delta_p2e + ...
        (transpose(grad_Delta_p2e)*grad_Delta_p2e).*del2_E2_del_Delta2;
    hess_E2_e2e = del_E2_del_Delta .* hess_Delta_e2e + ...
        (transpose(grad_Delta_e2e)*grad_Delta_e2e).*del2_E2_del_Delta2;
    hess_E2_p2t = del_E2_del_Delta .* hess_Delta_p2t + ...
        (transpose(grad_Delta_p2t)*grad_Delta_p2t).*del2_E2_del_Delta2;

    % generate matlab functions for gradient of E
    input = [transpose(x), Delta, h, K1];
    grad_E_shell_pen_p2p = matlabFunction(grad_E1_p2p, 'Vars', {input}, "File","grad_E_shell_pen_p2p.m");
    grad_E_shell_pen_p2e = matlabFunction(grad_E1_p2e, 'Vars', {input}, "File","grad_E_shell_pen_p2e.m");
    grad_E_shell_pen_e2e = matlabFunction(grad_E1_e2e, 'Vars', {input}, "File","grad_E_shell_pen_e2e.m");
    grad_E_shell_pen_p2t = matlabFunction(grad_E1_p2t, 'Vars', {input}, "File","grad_E_shell_pen_p2t.m");

    grad_E_shell_con_p2p = matlabFunction(grad_E2_p2p, 'Vars', {input}, "File","grad_E_shell_con_p2p.m");
    grad_E_shell_con_p2e = matlabFunction(grad_E2_p2e, 'Vars', {input}, "File","grad_E_shell_con_p2e.m");
    grad_E_shell_con_e2e = matlabFunction(grad_E2_e2e, 'Vars', {input}, "File","grad_E_shell_con_e2e.m");
    grad_E_shell_con_p2t = matlabFunction(grad_E2_p2t, 'Vars', {input}, "File","grad_E_shell_con_p2t.m");

    hess_E_shell_pen_p2p = matlabFunction(hess_E1_p2p, 'Vars', {input}, "File","hess_E_shell_pen_p2p.m");
    hess_E_shell_pen_p2e = matlabFunction(hess_E1_p2e, 'Vars', {input}, "File","hess_E_shell_pen_p2e.m");
    hess_E_shell_pen_e2e = matlabFunction(hess_E1_e2e, 'Vars', {input}, "File","hess_E_shell_pen_e2e.m");
    hess_E_shell_pen_p2t = matlabFunction(hess_E1_p2t, 'Vars', {input}, "File","hess_E_shell_pen_p2t.m");

    hess_E_shell_con_p2p = matlabFunction(hess_E2_p2p, 'Vars', {input}, "File","hess_E_shell_con_p2p.m");
    hess_E_shell_con_p2e = matlabFunction(hess_E2_p2e, 'Vars', {input}, "File","hess_E_shell_con_p2e.m");
    hess_E_shell_con_e2e = matlabFunction(hess_E2_e2e, 'Vars', {input}, "File","hess_E_shell_con_e2e.m");
    hess_E_shell_con_p2t = matlabFunction(hess_E2_p2t, 'Vars', {input}, "File","hess_E_shell_con_p2t.m");




end
