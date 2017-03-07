% test conic-ADMM3c vs. DE
clear; clc;
load ./tsp_test/tsplib
load ./qap_test/qaplib


%% QAPLIB
key = 'chr20a';
data = qaplib(key);
A = data(:,:,1); n = size(A,1);
B = data(:,:,2);

% n = 20;
% p = 0.1;
% noise = 0.05;
% A = binornd(1,p,n,n);
% B = -A-binornd(1,noise,n,n);

% define problem
problem = struct;
problem.A = A;          % graph 1
problem.B = B;          % graph 2
problem.P0 = eye(n);    % initialization
problem.maxits = 2000;  % maximum number of iterations
problem.varsize = 2;    % number of graph nodes per variable
problem.nprocs = 20;    % number of parallel processors

results = assign_csdp(problem);
results_de = assign_csdp_de(problem);
results_de132 = assign_csdp_de132(problem);

etas = results.eta_vec;
etas_de = results_de.eta_vec;
etas_de132 = results_de132.eta_vec;

plot(1:problem.maxits,log10(etas),1:problem.maxits,log10(etas_de),1:problem.maxits,log10(etas_de132),'LineWidth',1.5);
legend('Conic-ADMM3c (Sun et al.)','Direct Extension (1-2-3)', 'Direct Extension (1-3-2)', 'FontSize', 16);
ylabel('log10(\eta)', 'FontSize', 13)
xlabel('iterations', 'FontSize', 13)
set(gca, 'XTick', [0, 500, 1000, 1500, 2000])
set(gca,'FontSize',13)
print_figure('conic_vs_de_chr20a.pdf')

%% TSPLIB
key = 'gr21';
A = tsplib(key);
n = size(A,1);
B = diag(ones(1,n-1),1); B(1,n) = 1; B = 0.5*(B+B');

% define problem
problem = struct;
problem.A = A;          % graph 1
problem.B = B;          % graph 2
problem.P0 = eye(n);    % initialization
problem.maxits = 2000;  % maximum number of iterations
problem.varsize = 2;    % number of graph nodes per variable
problem.nprocs = 20;    % number of parallel processors

results = assign_csdp(problem);
results_de = assign_csdp_de(problem);
results_de132 = assign_csdp_de132(problem);

etas = results.eta_vec;
etas_de = results_de.eta_vec;
etas_de132 = results_de132.eta_vec;

plot(1:problem.maxits,log10(etas),1:problem.maxits,log10(etas_de),1:problem.maxits,log10(etas_de132),'LineWidth',1.5);
legend('Conic-ADMM3c (Sun et al.)','Direct Extension (1-2-3)', 'Direct Extension (1-3-2)', 'FontSize', 13);
ylabel('log10(\eta)', 'FontSize', 13)
xlabel('iterations', 'FontSize', 13)
set(gca, 'XTick', [0, 500, 1000, 1500, 2000])
set(gca,'FontSize',13)
print_figure('conic_vs_de_gr21.pdf')