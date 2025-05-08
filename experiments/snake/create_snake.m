% create snake like rod (sine wave)
N = 10;
xs = linspace(0, 2*pi, N);
ys = -sin(xs);

figure(1)
plot(xs, ys)

Nodes = [xs', ys', zeros(N,1)];
Edges = [(1:N-1)', (2:N)'];

filename = 'input_snake2.txt';
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