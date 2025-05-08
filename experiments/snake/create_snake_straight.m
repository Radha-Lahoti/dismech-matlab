% create straight rod of length 2pi
% Parameters
N = 10;                 % number of points
x_start = 0;             
x_end = 2*pi;            % range

x = linspace(x_start, x_end, N)';
y = zeros(N,1);
z = y;

Nodes = [x, y, z];
Edges = [(1:N-1)', (2:N)'];

filename = 'input_snake_straight.txt';
fid = fopen(filename,'w');
if fid ~= -1
    fprintf(fid,'*Nodes');
    fclose(fid);
end
writematrix(Nodes,filename, 'WriteMode','append')
fid = fopen(filename, 'at');
if fid ~= -1
    fprintf(fid,'*Edges');
    fclose(fid);
end
writematrix(Edges,filename, 'WriteMode','append')