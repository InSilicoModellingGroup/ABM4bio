clear, clc, close all

%cd /Users/vasvav/Public/invitro_neuro/examples/HeLa_cells/results/
data = table2array(readtable('results/stats.csv'));

T = data(:,1) / 24.0;
N_cells_pheno_0 = data(:, 4);
N_cells_pheno_1 = data(:, 5);
N_cells_pheno_1_Ap = data(:, 6);
N_cells_pheno_1_G1 = data(:, 7);
N_cells_pheno_1_Sy = data(:, 8);
N_cells_pheno_1_G2 = data(:, 9);
N_cells_pheno_1_Di = data(:,10);

N = N_cells_pheno_1;
G0_1 = N_cells_pheno_1_G1;
Sy   = N_cells_pheno_1_Sy;
G2_M = N_cells_pheno_1_G2 + N_cells_pheno_1_Di;
Ap = N_cells_pheno_1_Ap;

figure, hold on
plot(T,N), grid on
xlabel('days'), ylabel('number of cells')
hold off

figure, hold on
% plot(T, N, T, G0_1, T, Sy, T, G2_M), 
% legend('N', 'N: G0/1', 'N: S', 'N: G2/M')
plot(T, 100*G0_1./N, T, 100*Sy./N, T, 100*G2_M./N, T, 100*Ap./N), grid on 
xlabel('days'), ylabel('percentage of cells')
legend('%: G0/1', '%: S', '%: G2/M', '%: Ap')
hold off