function [Fc, Jc, Ffr, Jfr, shell_imc] = ...
    IMC_shell_force_only_trial(shell_imc, q, q0, dt)

k_c = shell_imc.k_c;
C = shell_imc.C;
delta = shell_imc.delta;
scale = shell_imc.scale;
% velTol = shell_imc.velTol;
contact_len = shell_imc.contact_len;
% mu_k = shell_imc.mu_k;
% friction_present = shell_imc.compute_friction;
n_dof = size(q,1);


if(~isempty(C))
    [colliding_tri_combos, colliding_constraint_types, ~,~] = detectShellCollisions_trial(q, C, constraint_type, delta, contact_len);

    use_hess = false; % if iter<omega compute only forces
    [Fc, Jc] = computeShellContact_trial...
        (q, colliding_tri_combos, colliding_constraint_types, delta, contact_len, scale, k_c, n_dof, use_hess);

    
    Ffr = zeros(n_dof, 1);
else 
    Fc = zeros(n_dof, 1);
    Ffr = zeros(n_dof, 1);
end

end
