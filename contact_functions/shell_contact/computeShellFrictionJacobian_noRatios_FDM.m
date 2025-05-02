function [Jfr, friction_partial_dfr_dfc, friction_partial_dfr_dx] = computeShellFrictionJacobian_noRatios_FDM(data, contact_forces, contact_jacobian, mu, dt, vel_tol, friction_type, constraint_type, ffr_original)

    change = 1e-8;

    friction_partial_dfr_dx = zeros(18,18);
    friction_partial_dfr_dfc = zeros(18,18);

    for i = 1:18
        changed_data = data;
        changed_data(i) = data(i)+change;

        [~, ffr_change] = computeShellFriction_noRatios(changed_data, contact_forces, mu, dt, vel_tol);
        friction_partial_dfr_dx (i,:) = (ffr_change - ffr_original) .* (1/change);
    end

    for j=1:18
        changed_forces = contact_forces;
        changed_forces(j) = contact_forces(j)+change;
        [~, ffr_change] = computeShellFriction_noRatios(data, changed_forces, mu, dt, vel_tol);
        friction_partial_dfr_dfc (j,:) = (ffr_change - ffr_original) .* (1/change);
    end

    Jfr = friction_partial_dfr_dfc*contact_jacobian + friction_partial_dfr_dx;
end

