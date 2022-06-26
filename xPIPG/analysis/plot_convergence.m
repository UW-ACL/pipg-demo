%{
05/18/2022
Purnanand Elango

Plot dual variable error for each exPIPG iteration until convergence for
feasible and infeasible problems
%}

clc
clear variables
close all

addpath ../
addpath ../solvers
addpath ../solvers/expipg
addpath ../utils

n = 64;
m = n/2;
N = 20;

x_max = 1;
u_max = 0.5;

pp = problem_data(n,m,N,x_max,u_max,1);

% Feasible problem
pp.expipg_feas_tol = 1e-7;

% Infeasible problem
% pp.expipg_infeas_tol = 1e-8;
% pp.z0(1:pp.n/2) = pp.z0(1:pp.n/2) + 0.7;



ppv = construct_vectorized(pp);
[~,~,~,err_dual,iter_conv] =  expipg_vec_diagnostic(ppv.xi{1},ppv.eta{1},...                                
                                                    ppv.xi{2},ppv.eta{2},ppv.xi{3},ppv.eta{3},100.1,100.1,100.1,...                
                                                    ppv.max_eig_P,pp.expipg_omg,pp.expipg_rho,pp.expipg_test_iter,pp.expipg_max_iter,-1,...  
                                                    pp.expipg_feas_tol,pp.expipg_infeas_tol,...                               
                                                    ppv.P,ppv.H,ppv.HT,...                                
                                                    ppv.xi_min,ppv.xi_max);


figure('Position',[100,100,400,400])
semilogy(1:iter_conv,err_dual(1:iter_conv),'-b','LineWidth',2);
title("$\|w^{j+1}-w^{j}\|_2$");
ylim([1e-3,1e3]);
xlim([1,0.65*iter_conv]);
xlabel("Iteration $j$");


rmpath ../
rmpath ../solvers
rmpath ../solvers/expipg
rmpath ../utils