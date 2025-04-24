function [grad_E1_p2p, grad_E1_p2e, grad_E1_e2e, grad_E2_p2p, grad_E2_p2e, grad_E2_e2e, ...
    hess_E1_p2p, hess_E1_p2e, hess_E1_e2e, hess_E2_p2p, hess_E2_p2e, hess_E2_e2e]...
    = generateContactHessianFunctions()

    % Declare symbolic variables
    syms x1s [1 3] real
    syms x1e [1 3] real
    syms x2s [1 3] real
    syms x2e [1 3] real
    syms Delta real
    syms h K1 real

    % point to point
    Delta_p2p = norm(x2s-x1s); % x1s and x2s are points of edge 1 and 2 in contact
    % we need to make sure that x1s and x2s are correctly put during simulation

    % point to edge
    e2 = x2e - x2s; % the edge for which the minimum distance vector does not lie on
    % an end
    e_other = x2s-x1s; % x1s is the node of edge 1 which the minimum distance vector does lie on
    Delta_p2e = norm(cross(e2,e_other))/norm(e2);

    % edge to edge
    e1 = x1e - x1s;
    u = cross(e1,e2);
    u_hat = u/norm(u);
    Delta_e2e = abs(dot((x1s-x2s),u_hat));

    % find gradient of Delta wrt [x_i, x_i+1, x_j, x_j+1]
    x = [x1s, x1e, x2s, x2e];
%     grad_Delta_p2p = jacobian(Delta_p2p, x);
grad_Delta_p2p = [x1s./Delta_p2p , zeros(1,3), x2s./Delta_p2p, zeros(1,3)];

    grad_Delta_p2e = jacobian(Delta_p2e, x);
    grad_Delta_e2e = jacobian(Delta_e2e, x);

    % find hessian of Delta wrt [x_i, x_i+1, x_j, x_j+1]
    hess_Delta_p2p = jacobian(grad_Delta_p2p, x);
    hess_Delta_p2e = jacobian(grad_Delta_p2e, x);
    hess_Delta_e2e = jacobian(grad_Delta_e2e, x);

    % Define Energy 
    % if Δ <= 2h - δ
    E1 = (2*h-Delta)^2;
    % if (2h - δ) < Δ < (2h + δ)
    E2 = ((1/K1)*log(1+exp(K1*(2*h-Delta))))^2;

    % find derivative of E wrt Delta
    del_E1_del_Delta = diff(E1,Delta);
    del_E2_del_Delta = diff(E2,Delta);

    % find second derivative of E wrt Delta
    del2_E1_del_Delta2 = diff(E1,Delta,2);
    del2_E2_del_Delta2 = diff(E2,Delta,2);

    % gradient of E
    grad_E1_p2p = del_E1_del_Delta .* grad_Delta_p2p;
    grad_E1_p2e = del_E1_del_Delta .* grad_Delta_p2e;
    grad_E1_e2e = del_E1_del_Delta .* grad_Delta_e2e;

    grad_E2_p2p = del_E2_del_Delta .* grad_Delta_p2p;
    grad_E2_p2e = del_E2_del_Delta .* grad_Delta_p2e;
    grad_E2_e2e = del_E2_del_Delta .* grad_Delta_e2e;

    % hessian of E
    hess_E1_p2p = del_E1_del_Delta .* hess_Delta_p2p + ...
        (transpose(grad_Delta_p2p)*grad_Delta_p2p).*del2_E1_del_Delta2;
    hess_E1_p2e = del_E1_del_Delta .* hess_Delta_p2e + ...
        (transpose(grad_Delta_p2e)*grad_Delta_p2e).*del2_E1_del_Delta2;
    hess_E1_e2e = del_E1_del_Delta .* hess_Delta_e2e + ...
        (transpose(grad_Delta_e2e)*grad_Delta_e2e).*del2_E1_del_Delta2;

    hess_E2_p2p = del_E2_del_Delta .* hess_Delta_p2p + ...
        (transpose(grad_Delta_p2p)*grad_Delta_p2p).*del2_E2_del_Delta2;
    hess_E2_p2e = del_E2_del_Delta .* hess_Delta_p2e + ...
        (transpose(grad_Delta_p2e)*grad_Delta_p2e).*del2_E2_del_Delta2;
    hess_E2_e2e = del_E2_del_Delta .* hess_Delta_e2e + ...
        (transpose(grad_Delta_e2e)*grad_Delta_e2e).*del2_E2_del_Delta2;


    % generate matlab functions for gradient of E
    input = [x, Delta, h, K1];
    grad_E_pen_p2p = matlabFunction(grad_E1_p2p, 'Vars', {input}, "File","grad_E_pen_p2p.m");
%     grad_E_pen_p2e = matlabFunction(grad_E1_p2e, 'Vars', {input}, "File","grad_E_pen_p2e.m");
%     grad_E_pen_e2e = matlabFunction(grad_E1_e2e, 'Vars', {input}, "File","grad_E_pen_e2e.m");

    grad_E_con_p2p = matlabFunction(grad_E2_p2p, 'Vars', {input}, "File","grad_E_con_p2p.m");
%     grad_E_con_p2e = matlabFunction(grad_E2_p2e, 'Vars', {input}, "File","grad_E_con_p2e.m");
%     grad_E_con_e2e = matlabFunction(grad_E2_e2e, 'Vars', {input}, "File","grad_E_con_e2e.m");

    hess_E_pen_p2p = matlabFunction(hess_E1_p2p, 'Vars', {input}, "File","hess_E_pen_p2p.m");
%     hess_E_pen_p2e = matlabFunction(hess_E1_p2e, 'Vars', {input}, "File","hess_E_pen_p2e.m");
%     hess_E_pen_e2e = matlabFunction(hess_E1_e2e, 'Vars', {input}, "File","hess_E_pen_e2e.m");

    hess_E_con_p2p = matlabFunction(hess_E2_p2p, 'Vars', {input}, "File","hess_E_con_p2p.m");
%     hess_E_con_p2e = matlabFunction(hess_E2_p2e, 'Vars', {input}, "File","hess_E_con_p2e.m");
%     hess_E_con_e2e = matlabFunction(hess_E2_e2e, 'Vars', {input}, "File","hess_E_con_e2e.m");


end
