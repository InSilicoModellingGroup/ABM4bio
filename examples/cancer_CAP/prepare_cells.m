%%% MATLAB code to prepare the cell data for the ABM4bio simulator
clc, clear

fid = fopen('cells.msh','r');
data = fscanf(fid,'%f %f %f %f',[4 Inf]);
fclose(fid);

nc = size(data,2);

fid = fopen('cells.dat','w');
fprintf(fid, '%d\n', nc);
for c = 1 : nc
    fprintf(fid, '%f %f %f\n', data(2,c),data(3,c),data(4,c));
    fprintf(fid, '%f\n', 15.0);
    %fprintf(fid, '%f\n', (25.0-15.0)*rand+15.0);
    fprintf(fid, '1.0 0.0 0.0\n0.0 1.0 0.0\n0.0 0.0 1.0\n');
    fprintf(fid, ' 1 0\n');
end
fclose(fid);
