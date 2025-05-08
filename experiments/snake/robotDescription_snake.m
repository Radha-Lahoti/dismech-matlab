% input: robotDescription.m

sim_params.static_sim = false;
sim_params.TwoDsim = true;
sim_params.use_midedge = false; % boolean var to decide on using midedge normal or 
% hinge model for shell bending
sim_params.use_lineSearch = false;
sim_params.log_data = true;
sim_params.logStep = 10;
sim_params.showFrames = false;

% Time step
sim_params.dt = 1e-2;

% Maximum number of iterations in Newton Solver
sim_params.maximum_iter = 100;

% Total simulation time
if(sim_params.static_sim)
%     sim_params.totalTime = sim_params.dt;
    sim_params.totalTime = sim_params.dt*10;
else
    sim_params.totalTime = 15; % sec
end

% How often the plot should be saved? (Set plotStep to 1 to show each plot)
sim_params.plotStep = 100;

%% Input text file 
% inputFileName = 'experiments/snake/input_snake1.txt';
% inputFileName = 'experiments/snake/input_snake2.txt';
inputFileName = 'experiments/snake/input_snake_straight.txt';

% reading the input text file
[nodes, edges, face_nodes] = inputProcessorNew(inputFileName);
%% Input parameters
% geometry parameters
geom.rod_r0 = 1e-3;
geom.shell_h = 1e-3;
% 
% material parameters
material.density = 1500;
material.youngs_rod = 2e9;
material.youngs_shell = 0;
material.poisson_rod = 0.3;
material.poisson_shell = 0.5;
% 
%% external force list ["selfContact", "selfFriction", "floorContact", "floorFriction", "gravity", "buoyancy", "viscous", "aerodynamic","pointForce"]
env.ext_force_list = ["gravity", "viscous"]; 

% environment parameters
env.g = [0, 0, 0]';
env.eta = 0.0;

%% Tolerance on force function. 

sim_params.tol = 1e-4;
sim_params.ftol = 1e-4;
sim_params.dtol = 1e-2;

%% Boundary conditions
fixed_node_indices = [];
fixed_edge_indices = [];
input_log_node = 1;

%% initial conditions

%% Plot dimensions
sim_params.plot_x = [-0.1,2*pi];
sim_params.plot_y = [-2,2];
sim_params.plot_z = [-0.1,0.1];
sim_params.view = [0,90]; % x-y plane
