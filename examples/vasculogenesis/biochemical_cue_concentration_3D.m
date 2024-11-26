clear, clc, close all

S_min = -100.0;
S_max = -S_min;
N = 50;
DS = (S_max - S_min) / N;

R_max = sqrt(3.0)*(S_max-S_min)/2.0;

Field_min = 0.020e-0;
Field_max = 0.980e-0;
noise_amp = 0.020e-0;

fileID = fopen('VEGF.radial.dat', 'w');

fprintf(fileID, '%d\n',(N-1)^3);
for I = 1 : N-1
    X = S_min + DS * I;
    for J = 1 : N-1
        Y = S_min + DS * J;
        for K = 1 : N-1
            Z = S_min + DS * K;
            R = sqrt(X.^2 + Y.^2 + Z.^2);
            %
            a = -(Field_max-Field_min);
            b = Field_max;
            Field = a * (R/R_max) + b;
            % add some noise
            Field_noise = noise_amp * (-1.0 + 2.0*rand);
            Field = Field + Field_noise;
            fprintf(fileID, '%f %f %f %f\n', X,Y,Z,Field);
        end
    end
end

fclose(fileID);
