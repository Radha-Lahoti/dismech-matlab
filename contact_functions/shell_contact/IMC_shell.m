function [Fc, Jc, Ffr, Jfr, shell_imc] = ...
    IMC_shell(shell_imc, q, q0, tri_combos, iter, dt, fvec_exceptIMC, fixedDOF)

k_c = shell_imc.k_c;
C = shell_imc.C;
delta = shell_imc.delta;
omega = shell_imc.omega;
scale = shell_imc.scale;
velTol = shell_imc.velTol;
contact_len = shell_imc.contact_len;
mu_k = shell_imc.mu_k;
friction_present = shell_imc.compute_friction;


col_lim = 100*delta;
candidate_lim = scale*(contact_len + col_lim);
n_dof = size(q,1);

if(iter==1) % run only on first iter
    [C, ~] = constructShellCandidateSet(q, tri_combos, candidate_lim, scale);
    if(~isempty(C)) % if collision is detected, update contact stiffness if necessary
        k_c = updateContactStiffnessNew(fvec_exceptIMC, C, fixedDOF);
    end
end
if(~isempty(C))
    [colliding_tri_combos, ~,~] = detectShellCollisions(q, C, delta, contact_len);

    use_hess = false; % if iter<omega compute only forces
    if (iter>omega) 
        use_hess = true; % compute Jacobian for convergence
    end
%     [Fc, Jc] = computeShellContact...
%         (q, colliding_tri_combos, delta, contact_len, scale, k_c, n_dof, use_hess);

%     Ffr = zeros(n_dof, 1);
%     Jfr = zeros(n_dof, n_dof);

    [Fc, Jc, Ffr, Jfr] = computeShellContactAndFriction...
        (q, q0, colliding_tri_combos, delta, contact_len, scale, k_c, mu_k, dt, velTol, n_dof, use_hess, friction_present);


else 
    Fc = zeros(n_dof, 1);
    Jc = zeros(n_dof, n_dof);
    Ffr = zeros(n_dof, 1);
    Jfr = zeros(n_dof, n_dof);
end
shell_imc.C= C;
shell_imc.k_c = k_c;

end
