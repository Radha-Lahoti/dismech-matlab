function [friction_partial_dfr_dx1_func, friction_partial_dfr_dfc1_func, ...
    friction_partial_dfr_dx2_func, friction_partial_dfr_dfc2_func] ...
    = generate_shell_friction_jacobian_function()
    % Generate Jacobian of friction force.
    % Treat friction direction explicitly using previous time step's velocity.
    % Contact force magnitude is treated implicitly using contact energy Hessian.

    disp('Starting friction Jacobian...')

    % Define symbolic variables
    syms mu dt vel_tol real
    syms x11 [3 1] real
    syms x12 [3 1] real
    syms x13 [3 1] real
    syms x21 [3 1] real
    syms x22 [3 1] real
    syms x23 [3 1] real
    syms x11_0 [3 1] real
    syms x12_0 [3 1] real
    syms x13_0 [3 1] real
    syms x21_0 [3 1] real
    syms x22_0 [3 1] real
    syms x23_0 [3 1] real

    syms fc11 [3 1] real
    syms fc12 [3 1] real
    syms fc13 [3 1] real
    syms fc21 [3 1] real
    syms fc22 [3 1] real
    syms fc23 [3 1] real

    syms ratios [6 1] real

    K2 = 15/vel_tol;

    % Compute norms
    fc11n = norm(fc11);
    fc12n = norm(fc12);
    fc13n = norm(fc13);
    fc21n = norm(fc21);
    fc22n = norm(fc22);
    fc23n = norm(fc23);

    % ratios
    r1 = ratios(1); s1 = ratios(2); t1 = ratios(3);
    r2 = ratios(4); s2 = ratios(5); t2 = ratios(6);
%% Compute relative velocities
    % tri1
    v11 = (x11 - x11_0)./dt;
    v12 = (x12 - x12_0)./dt;
    v13 = (x13 - x13_0)./dt;
    % tri2
    v21 = (x21 - x21_0)./dt;
    v22 = (x22 - x22_0)./dt;
    v23 = (x23 - x23_0)./dt;

    v1 = r1*v11 + s1*v12 + t1*v13;
    v2 = r2*v21 + s2*v22 + t2*v23;
    v_rel = v1 - v2;

    % direction of contact
    p1 = x11*r1+x12*s1+x13*t1; % point of contact on triangle 1
    p2 = x21*r2+x22*s2+x23*t2; % point of contact on triangle 2
    contact_dir = p1 - p2;
%     assert(norm(contact_dir)>0, "contact distance is 0, can't compute direction of contact");
    contact_norm = contact_dir/norm(contact_dir);

    % Tangential relative velocity
    tv_rel = v_rel - (dot(v_rel, contact_norm).* contact_norm);
    tv_rel_n = norm(tv_rel);
    tv_rel_u = tv_rel / tv_rel_n;

    gamma = (2 / (1 + exp(-K2*tv_rel_n))) - 1;

    % STICKING FRICTION JACOBIAN 0 < γ < 1
    ffr1 = (gamma * mu) .* tv_rel_u;
    ffr2 = -ffr1;

    ffr11 = ffr1 .* fc11n;
    ffr12 = ffr1 .* fc12n;
    ffr13 = ffr1 .* fc13n;
    ffr21 = ffr2 .* fc21n;
    ffr22 = ffr2 .* fc22n;
    ffr23 = ffr2 .* fc23n;

    ffr_vec1 = [ffr11; ffr12; ffr13; ffr21; ffr22; ffr23];

    wrt_nodes = [x11; x12; x13; x21; x22; x23];
    wrt_cforces = [fc11; fc12; fc13; fc21; fc22; fc23];

    friction_partial_dfr_dx1 = jacobian(ffr_vec1, wrt_nodes);
    friction_partial_dfr_dfc1 = jacobian(ffr_vec1, wrt_cforces);

    % SLIDING FRICTION JACOBIAN γ >= 1
    ffr1 = mu .* tv_rel_u;
    ffr2 = -ffr1;

    ffr11 = ffr1 .* fc11n;
    ffr12 = ffr1 .* fc12n;
    ffr13 = ffr1 .* fc13n;
    ffr21 = ffr2 .* fc21n;
    ffr22 = ffr2 .* fc22n;
    ffr23 = ffr2 .* fc23n;

    ffr_vec2 = [ffr11; ffr12; ffr13; ffr21; ffr22; ffr23];

    friction_partial_dfr_dx2 = jacobian(ffr_vec2, wrt_nodes);
    friction_partial_dfr_dfc2 = jacobian(ffr_vec2, wrt_cforces);

    % Combine inputs for matlab function generation
    inputs = [x11; x12; x13; x21; x22; x23; x11_0; x12_0; x13_0; x21_0; x22_0; x23_0; fc11; fc12; fc13; fc21; fc22; fc23; ratios; mu; dt; vel_tol];

    % Convert to different MATLAB functions
%     friction_partial_dfr_dx1_func = matlabFunction(friction_partial_dfr_dx1, 'Vars', {inputs}, "File","stickFriction_partial_dfr_dx_func.m");
%     friction_partial_dfr_dfc1_func = matlabFunction(friction_partial_dfr_dfc1, 'Vars', {inputs}, "File","stickFriction_partial_dfr_dfc_func.m");
%     friction_partial_dfr_dx2_func = matlabFunction(friction_partial_dfr_dx2, 'Vars', {inputs}, "File","slideFriction_partial_dfr_dx_func.m");
%     friction_partial_dfr_dfc2_func = matlabFunction(friction_partial_dfr_dfc2, 'Vars', {inputs}, "File","slideFriction_partial_dfr_dfc_func.m");

    disp('Completed friction Jacobian calculation.');

end
