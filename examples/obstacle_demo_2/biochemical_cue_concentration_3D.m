clear, clc, close all

S_min = -100.0;
S_max = -S_min;
N = 50;
DS = (S_max - S_min) / N;

R_max = sqrt(3.0)*(S_max-S_min)/2.0;

Field_min = 0.999e-0;
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
            a = (Field_min-Field_max);
            b = Field_max;
            Field = a * (R/R_max) + b;
            fprintf(fileID, '%f %f %f %f\n', X,Y,Z,Field);
        end
    end
end

fclose(fileID);

fileID = fopen('O2.axial.dat', 'w');

fprintf(fileID, '%d\n',(N-1)^3);
for I = 1 : N-1
    X = S_min + DS * I;
    for J = 1 : N-1
        Y = S_min + DS * J;
        for K = 1 : N-1
            Z = S_min + DS * K;
            %
            a = -(Field_max-Field_min)/(S_max-S_min);
            b = Field_min - a*S_min;
            Field = a * X + b;
            fprintf(fileID, '%f %f %f %f\n', X,Y,Z,Field);
        end
    end
end

fclose(fileID);

fileID = fopen('O2.radial2axial.dat', 'w');

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
            if ((S_max+S_min)/2.0) < X
                a = (Field_min-Field_max);
                b = Field_max;
                Field = a * (R/R_max) + b;
                fprintf(fileID, '%f %f %f %f\n', X,Y,Z,Field);
            else
                a = -(Field_max-Field_min)/(S_max-S_min);
                b = Field_min - a*S_min;
                Field = a * X + b;
                fprintf(fileID, '%f %f %f %f\n', X,Y,Z,Field);
            end
        end
    end
end

fclose(fileID);
