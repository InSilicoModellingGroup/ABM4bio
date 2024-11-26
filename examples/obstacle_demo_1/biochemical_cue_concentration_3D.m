clear, clc, close all

S_min = -100.0;
S_max = -S_min;
N = 50;
DS = (S_max - S_min) / N;

R_max = sqrt(3.0)*(S_max-S_min)/2.0;

Field_min = 0.000e-0;
Field_max = 1.000e-0;

fileID = fopen('O2.radial.dat', 'w');

fprintf(fileID, '%d\n',(N-1)^3);
for I = 1 : N-1
    X = S_min + DS * I;
    for J = 1 : N-1
        Y = S_min + DS * J;
        for K = 1 : N-1
            Z = S_min + DS * K;
            %
            R = sqrt(X.^2 + Y.^2 + Z.^2);
            %
            a = (Field_max-Field_min);
            b = Field_min;
            Field = a * (R/R_max) + b;
            fprintf(fileID, '%f %f %f %f\n', X,Y,Z,Field);
        end
    end
end

fclose(fileID);

T = 100;

Field_max = 1.00e-0;

system('rm -Rf Velocity/; mkdir Velocity');
for time = 0 : T
    
    file_name = sprintf("Velocity/%d.dat", time);
    fileID = fopen(file_name, 'w');
    
    fprintf(fileID, '%d\n',(N-1)^3);
    for I = 1 : N-1
        X = S_min + DS * I;
        for J = 1 : N-1
            Y = S_min + DS * J;
            for K = 1 : N-1
                Z = S_min + DS * K;
                R = sqrt(X.^2 + Y.^2 + Z.^2);
                %
                Field = (Field_max/R) * cos(time*2.0*pi/T);
                fprintf(fileID, '%f %f %f %f %f %f\n', X,Y,Z,Field*X,Field*Y,Field*Z);
            end
        end
    end

    fclose(fileID);

end
