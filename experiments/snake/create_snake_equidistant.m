clc 
clear all
% Parameters
N = 10;                 % number of points
x_start = 0;             
x_end = 2*pi;            % range over which to draw the sine curve

% High-resolution x for arc length computation
x_fine = linspace(x_start, x_end, 10000);
y_fine = sin(x_fine);

% Approximate arc length
dx = diff(x_fine);
dy = diff(y_fine);
ds = sqrt(dx.^2 + dy.^2);          % segment lengths
s = [0, cumsum(ds)];               % cumulative arc length

% Desired equally spaced arc lengths
s_equidistant = linspace(0, s(end), N);

% Interpolate to get x and y at equally spaced arc lengths
x_equidistant = interp1(s, x_fine, s_equidistant);
y_equidistant = -sin(x_equidistant);

% Store into nodes array
nodes = [x_equidistant(:), y_equidistant(:)];

% Optional: visualize
plot(x_fine, y_fine, 'k--'); hold on;
plot(nodes(:,1), nodes(:,2), 'ro-'); axis equal;

tolerance = 1e-10;
nodes(abs(nodes) < tolerance) = 0;
nodes = [nodes, zeros(N,1)];
Edges = [(1:N-1)', (2:N)'];

%% put into inputGeometry.txt file
filename = 'input_snake2.txt';
fid = fopen(filename,'w');
if fid ~= -1
    fprintf(fid,'*Nodes');
    fclose(fid);
end
writematrix(nodes,filename, 'WriteMode','append')
fid = fopen(filename, 'at');
if fid ~= -1
    fprintf(fid,'*Edges');
    fclose(fid);
end
writematrix(Edges,filename, 'WriteMode','append')
