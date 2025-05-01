% input: robotDescription.m

sim_params.static_sim = false;
sim_params.TwoDsim = false;
sim_params.use_midedge = false; % boolean var to decide on using midedge normal or 
% hinge model for shell bending
sim_params.use_lineSearch = 0;
sim_params.showFrames = false;
sim_params.logStep = 1;
sim_params.log_data = true;

% Time step
sim_params.dt = 1e-3;

% Maximum number of iterations in Newton Solver
sim_params.maximum_iter = 100;

% Total simulation time
if(sim_params.static_sim)
%     sim_params.totalTime = sim_params.dt;
    sim_params.totalTime = sim_params.dt*5;
else
    sim_params.totalTime = 2; % sec
end

% How often the plot should be shown? (Set plotStep to 1 to show each plot)
sim_params.plotStep = 10;

%% Input parameters
% geometry parameters
geom.rod_r0 = 0;
geom.shell_h = 5e-3;

% material parameters
material.density = 1200;
material.youngs_rod = 0; % not used
material.youngs_shell = 2e8;
material.poisson_rod = 0;
material.poisson_shell = 0.5;

%% external force list ["selfContact", "selfFriction", "floorContact", "floorFriction", "gravity", "buoyancy", "viscous", "aerodynamic","pointForce"]
env.ext_force_list = ["selfContact", "selfFriction"]; 
% env.ext_force_list = ["selfContact"]; 

% environment parameters
env.g = [0, 0, -9.81]';
material.contact_stiffness = 0.1;
env.velTol = 1e-2;
material.mu = 0.25;

%% Input text file 
inputFileName = 'experiments/shellContact/input_twoTriangleContact.txt';
% inputFileName = 'experiments/shellContact/input_twoTriangleContact_p2p.txt';
% inputFileName = 'experiments/shellContact/input_twoTriangleContact_p2e.txt';
% inputFileName = 'experiments/shellContact/input_twoTriangleContact_p2t.txt';
% inputFileName = 'experiments/shellContact/input_twoTriangleContact_e2e.txt';

% reading the input text file
[nodes, edges, face_nodes] = inputProcessorNew(inputFileName);

%% Tolerance on force function. 

sim_params.tol = 1e-4;
sim_params.ftol = 1e-4;
sim_params.dtol = 1e-4;

%% Boundary conditions
fixed_node_indices = [];
fixed_edge_indices = [];

%% logging
input_log_node = 1;

%% initial conditions
u_init = -0.5;

%% Plot dimensions
sim_params.plot_x = [-2,2];
sim_params.plot_y = [-2,2];
sim_params.plot_z = [-2,2];
