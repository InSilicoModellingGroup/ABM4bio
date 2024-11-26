clear, clc, close all

S_min = -100.0;
S_max = -S_min;
N = 50;
DS = (S_max - S_min) / N;

fileID = fopen('field_Convection.dat', 'w');

fprintf(fileID, '%d\n',(N-1)^3);
for I = 1 : N-1
    X = S_min + DS * I;
    for J = 1 : N-1
        Y = S_min + DS * J;
        for K = 1 : N-1
            Z = S_min + DS * K;
            %
            R = sqrt(X.^2+Y.^2);
            %
            if R < 40.0
                Vz = 3.0;
            else
                Vz = 0.0;
            end
            fprintf(fileID, '%f %f %f 0. 0. %f\n', X,Y,Z,Vz);
        end
    end
end

fclose(fileID);

Field_min = 0.890e-0;
Field_max = 1.000e-0;

fileID = fopen('field_O2.dat', 'w');

fprintf(fileID, '%d\n',(N-1)^3);
for I = 1 : N-1
    X = S_min + DS * I;
    for J = 1 : N-1
        Y = S_min + DS * J;
        for K = 1 : N-1
            Z = S_min + DS * K;
            %
            R = sqrt(X.^2+Y.^2);
            %
            a = (Field_min-Field_max);
            b = Field_max;
            if R < 40.0
                Field = Field_max;
            elseif R < 50.0
                Field = a * ((R-40.0)/(50.0-40.0)) + b;
            else
                Field = Field_min;
            end
            fprintf(fileID, '%f %f %f %f\n', X,Y,Z,Field);
        end
    end
end

fclose(fileID);

Field_min = 0.900e-0;
Field_max = 1.000e-0;

fileID = fopen('field_ECM.dat', 'w');

fprintf(fileID, '%d\n',(N-1)^3);
for I = 1 : N-1
    X = S_min + DS * I;
    for J = 1 : N-1
        Y = S_min + DS * J;
        for K = 1 : N-1
            Z = S_min + DS * K;
            %
            R = sqrt(X.^2+Y.^2);
            %
            if R < 40.0
                Field = 0;
            else
                Field = Field_max;%Field_min+(Field_max-Field_min)*rand;
            end
            fprintf(fileID, '%f %f %f %f\n', X,Y,Z,Field);
        end
    end
end

fclose(fileID);
