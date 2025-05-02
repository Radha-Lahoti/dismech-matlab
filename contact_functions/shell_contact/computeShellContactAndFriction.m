function [Fc, Jc, Ffr, Jfr] = computeShellContactAndFriction(q, q0,  C, delta, contact_len, scale, k_c, mu_k, dt, vel_tol, n_dof, use_hess, friction_present)
% % Inputs
%   C: triangle_combos which can potentially come in contact: Candidate set
% % Outputs
%   Fc: contact force
% -----------------------------------
assert(size(C, 2) == 6);
num_inputs = size(C, 1);
Fc = zeros(n_dof,1);
Jc = zeros(n_dof,n_dof);
Ffr = zeros(n_dof,1);
Jfr = zeros(n_dof,n_dof);

K1 = 15*contact_len/(2*delta);
contact_lim   = scale*(contact_len + delta);
numerical_lim = scale*(contact_len - delta);

for i = 1:num_inputs
    [dist, constraint_type, tri_combo_input, combo_nodes_updated, ratios] = modified_triangleContactType(q, C(i,:), scale);
%     dist
%     constraint_type
%     numerical_lim
    %% Contact
    % input
    input = [tri_combo_input', dist, contact_len/2*scale, K1];

    if ( dist<=numerical_lim )
        % if Δ <= 2h - δ: penetration
        if (constraint_type=="PointToPoint")
            gradEc = grad_E_shell_pen_p2p(input);

            if(use_hess)
                hessEc = hess_E_shell_pen_p2p(input);
            else
                hessEc = zeros(18,18);
            end

        elseif (constraint_type=="PointToEdge")
            gradEc = grad_E_shell_pen_p2e(input);

            if(use_hess)
                hessEc = hess_E_shell_pen_p2e(input);
            else
                hessEc = zeros(18,18);
            end

        elseif (constraint_type=="EdgeToEdge")
            gradEc = grad_E_shell_pen_e2e(input);

            if(use_hess)
                hessEc = hess_E_shell_pen_e2e(input);
            else
                hessEc = zeros(18,18);
            end

        elseif (constraint_type=="PointToFace")
            gradEc = grad_E_shell_pen_p2t(input);

            if(use_hess)
                hessEc = hess_E_shell_pen_p2t(input);
            else
                hessEc = zeros(18,18);
            end
        elseif (constraint_type=="EdgeToFace")
            gradEc = zeros(1,18);
            hessEc = zeros(18,18);
        end

    elseif ( dist > numerical_lim && dist < contact_lim )
        % if (2h - δ) < Δ < (2h + δ): contact zone but no penetration
        if (constraint_type=="PointToPoint")
            gradEc = grad_E_shell_con_p2p(input);

            if(use_hess)
                hessEc = hess_E_shell_con_p2p(input);
            else
                hessEc = zeros(18,18);
            end

        elseif (constraint_type=="PointToEdge")
            gradEc = grad_E_shell_con_p2e(input);

            if(use_hess)
                hessEc = hess_E_shell_con_p2e(input);
            else
                hessEc = zeros(18,18);
            end

        elseif (constraint_type=="EdgeToEdge")
            gradEc = grad_E_shell_con_e2e(input);

            if(use_hess)
                hessEc = hess_E_shell_con_e2e(input);
            else
                hessEc = zeros(18,18);
            end
        elseif (constraint_type=="PointToFace")
            gradEc = grad_E_shell_con_p2t(input);

            if(use_hess)
                hessEc = hess_E_shell_con_p2t(input);
            else
                hessEc = zeros(18,18);
            end
        elseif (constraint_type=="EdgeToFace")
            gradEc = zeros(1,18);
            hessEc = zeros(18,18);
        end
    else
        gradEc = zeros(1,18);
        hessEc = zeros(18,18);
    end
    fc = (scale * k_c).*gradEc;
    jc = (scale^2 * k_c).*hessEc;

    ind = [mapNodetoDOF(combo_nodes_updated(1)); mapNodetoDOF(combo_nodes_updated(2)); mapNodetoDOF(combo_nodes_updated(3));...
        mapNodetoDOF(combo_nodes_updated(4)); mapNodetoDOF(combo_nodes_updated(5)); mapNodetoDOF(combo_nodes_updated(6))];

    Fc(ind) = Fc(ind) - fc';
    Jc(ind,ind) = Jc(ind,ind) - jc;

        %% Friction
        if(friction_present)
            for j = 1:6
                tri_combo_input0(3*j-2:3*j) = q0(ind (3*j-2:3*j));
            end
    
            data = [tri_combo_input', tri_combo_input0];
            if (find(gradEc)) % only compute friction force if there is non-zero contact force
%                 [friction_type,ffr] = computeShellFriction(data, ratios, fc, mu_k, dt, vel_tol);
                [friction_type,ffr] = computeShellFriction_noRatios(data, fc, mu_k, dt, vel_tol);

                Ffr(ind) = Ffr(ind) - ffr';
    
                if(use_hess)
                    if(friction_type=="ZeroVel")
                        jfr = zeros(18,18);
                    else
%                         jfr = computeShellFrictionJacobian(data, ratios, fc, jc, mu_k, dt, vel_tol, friction_type, constraint_type);
%                         jfr = computeShellFrictionJacobian_noRatios(data, fc, jc, mu_k, dt, vel_tol, friction_type, constraint_type);
                        jfr = computeShellFrictionJacobian_noRatios_FDM(data, fc, jc, mu_k, dt, vel_tol, friction_type, constraint_type, ffr);
    
                    end
                    Jfr(ind,ind) = Jfr(ind,ind) - jfr;
                end
            end
        end

end
end


