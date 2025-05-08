clc
clear all
Tp = 2*pi;
totalT = 15;
omega = 2*pi/Tp;

dt = 1e-2;

Nsteps = round(totalT/dt);
ctime = 0; % current time
time_arr = linspace(0,totalT,Nsteps);

filename1 = 'kappa_bar1.xls';
filename2 = 'kappa_bar2.xls';
kappa_bar1 = readmatrix(filename1);
kappa_bar2 = readmatrix(filename2);

amplitudes = kappa_bar2 - kappa_bar1;

n_springs = size(amplitudes,1);

% phase difference between front end and back end
total_phase_difference = 0; 
del_phase = total_phase_difference/(n_springs-1);


k_bar = zeros(n_springs, 2, Nsteps);

for timeStep = 1:Nsteps

    for i=1:n_springs
        phase = (i-1)*del_phase;
        k_bar(i,2,timeStep) = amplitudes(i)/2*(cos(omega*ctime + phase) - 1);
    end

    ctime = ctime + dt;
end

save('kappa_bar.mat', 'k_bar');