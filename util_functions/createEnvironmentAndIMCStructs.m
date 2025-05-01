function [environment,imc, shell_imc] = createEnvironmentAndIMCStructs(env,geom,material,sim_params)

environment = struct();
imc = struct();
shell_imc = struct();

environment.ext_force_list = env.ext_force_list;

if ismember("gravity",env.ext_force_list)
   environment.g = env.g;
end

if ismember("buoyancy",env.ext_force_list)
    environment.rho = env.rho;
end

if ismember("viscous", env.ext_force_list)
    environment.eta = env.eta;
end

if ismember("aerodynamic", env.ext_force_list)
    environment.rho = env.rho;
    environment.Cd = env.Cd;
end

if ismember("pointForce", env.ext_force_list)
    environment.ptForce = env.ptForce;
    environment.ptForce_node = env.ptForce_node;
end

if ismember("selfContact", env.ext_force_list)
    imc.k_c = material.contact_stiffness;
    imc.contact_len = 2*geom.rod_r0;
    imc.delta = 0.01*imc.contact_len;
    imc.omega = 20; % # iters before jacobian for contact forces is used
    imc.scale = 1/geom.rod_r0;
    imc.C = [];

    % shell contact
    shell_imc.k_c = material.contact_stiffness;
    shell_imc.contact_len = 2*geom.shell_h;
    shell_imc.delta = 0.01*shell_imc.contact_len;
    shell_imc.omega = 20; % # iters before jacobian for contact forces is used
    shell_imc.scale = 1/geom.shell_h;
    shell_imc.C = [];
    shell_imc.constraint_type = string;

    if ismember("selfFriction", env.ext_force_list)
        imc.compute_friction = true;
        imc.mu_k = material.mu;
        imc.velTol = env.velTol;

        % shell contact
        shell_imc.compute_friction = true;
        shell_imc.mu_k = material.mu;
        shell_imc.velTol = env.velTol;
    else
        imc.compute_friction = false;
        shell_imc.compute_friction = false;

        imc.mu_k = [];
        imc.velTol = [];
        shell_imc.mu_k = [];
        shell_imc.velTol = [];
    end
end

if ismember("floorContact", env.ext_force_list)
    imc.floor_z = env.floor_z;
    imc.k_c_floor = env.contact_stiffness;
    imc.h = geom.rod_r0;
    imc.delta_floor = 1*geom.rod_r0;
    imc.omega = 20; % # iters before jacobian for contact forces is used
    imc.scale = 1/imc.h;
    environment.showFloor = true;

    if ismember("floorFriction", env.ext_force_list)
        imc.floor_has_friction = true;
        imc.mu_floor = env.mu;
        imc.velTol = env.velTol;
    else
        imc.floor_has_friction = false;
    end
end
if sim_params.static_sim
    environment.static_g = env.g;
end
