function Jfr = computeShellFrictionJacobian_noRatios(data, contact_forces, contact_jacobian, mu, dt, vel_tol, friction_type, constraint_type)
    % Extract position data
    x11 = data(1:3);
    x12 = data(4:6);
    x13 = data(7:9);
    x21 = data(10:12);
    x22 = data(13:15);
    x23 = data(16:18);
    x11_0 = data(19:21);
    x12_0 = data(22:24);
    x13_0 = data(24:27);
    x21_0 = data(28:30);
    x22_0 = data(31:33);
    x23_0 = data(34:36);

    % Extract force data
    fc11 = contact_forces(1:3);
    fc12 = contact_forces(4:6);
    fc13 = contact_forces(7:9);
    fc21 = contact_forces(10:12);
    fc22 = contact_forces(13:15);
    fc23 = contact_forces(16:18);

    inputs = [x11, x12, x13, x21, x22, x23, x11_0, x12_0, x13_0, x21_0, x22_0, x23_0, fc11, fc12, fc13, fc21, fc22, fc23, mu, dt, vel_tol]';
    
    if(friction_type == "Sticking")
        friction_partial_dfr_dx = stickFriction_partial_dfr_dx_func(inputs);
        friction_partial_dfr_dfc = stickFriction_partial_dfr_dfc_func(inputs);
%         friction_partial_dfr_dx = slideFriction_partial_dfr_dx_func(inputs);
%         friction_partial_dfr_dfc = slideFriction_partial_dfr_dfc_func(inputs);

    elseif(friction_type == "Sliding")
%         friction_partial_dfr_dx = stickFriction_partial_dfr_dx_func(inputs);
%         friction_partial_dfr_dfc = stickFriction_partial_dfr_dfc_func(inputs);
        friction_partial_dfr_dx = slideFriction_partial_dfr_dx_func(inputs);
        friction_partial_dfr_dfc = slideFriction_partial_dfr_dfc_func(inputs);
    else
        error("friction_type for Jacobian computation should be either Sticking or Sliding");
    end

    if constraint_type == "PointToPoint"
        friction_partial_dfr_dfc(:, 4:9) = zeros(18,6);  % remove the entries corresponding to x12, x13
        friction_partial_dfr_dfc(:, 13:18) = zeros(18,6); % remove the entries corresponding to x22, x23
    elseif constraint_type == "PointToEdge"
        friction_partial_dfr_dfc(:, 4:9) = zeros(18,6);  % remove the entries corresponding to x12, x13
        friction_partial_dfr_dfc(:, 16:18) = zeros(18,3); % remove the entries corresponding to x23
    elseif constraint_type == "EdgeToEdge"
        friction_partial_dfr_dfc(:, 7:9) = zeros(18,3);  % remove the entries corresponding to x13
        friction_partial_dfr_dfc(:, 16:18) = zeros(18,3); % remove the entries corresponding to x23
    elseif constraint_type == "PointToTriangle"
        friction_partial_dfr_dfc(:, 4:9) = zeros(18,6);  % remove the entries corresponding to x12, x13
    elseif constraint_type == "EdgeToTriangle"
        friction_partial_dfr_dfc(:, 7:9) = zeros(18,3);  % remove the entries corresponding to x13
    end

    
    Jfr = friction_partial_dfr_dfc*contact_jacobian + friction_partial_dfr_dx;
%     assert(~anynan(Jfr),'IMC friction jacobian is not real (NaN).');


end