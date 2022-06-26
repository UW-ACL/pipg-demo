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

n = 32;
m = n/2;
N = 20;

x_max = 1;
u_max = 0.5;

omg = 200;
alf = nan;

pp = problem_data(n,m,N,x_max,u_max,1);
pp.expipg_omg = omg;
ppv = construct_vectorized(pp);
lam = ppv.max_eig_P;
sig = ppv.spec_norm_HTH;

rho_vec = [0.1, 1.0, 1.5, 1.9];

err = cell(length(rho_vec), 1);
iter_conv = zeros(length(rho_vec), 1);

for k = 1:length(rho_vec)

    alf = 2/((lam + 4*omg*sig)^0.5+lam);
    bet = omg*alf;

    [~,~,~,err{k},iter_conv(k)] =  expipg_vec_diagnostic(ppv.xi{1},ppv.eta{1},...
                                   ppv.xi{2},ppv.eta{2},ppv.xi{3},ppv.eta{3},100.1,100.1,100.1,...
                                   lam,omg,rho_vec(k),pp.expipg_test_iter,pp.expipg_max_iter,-1,...
                                   pp.expipg_feas_tol,pp.expipg_infeas_tol,...
                                   ppv.P,ppv.H,ppv.HT,...
                                   ppv.xi_min,ppv.xi_max);
    err{k} = err{k}./(rho_vec(k)*bet);
    
end

%%

min_iter_conv = min(iter_conv);
figure('Position',[100,100,400,400])
semilogy(1:min_iter_conv,err{1}(1:min_iter_conv),'-g','DisplayName',num2str(rho_vec(1)));
hold on
semilogy(1:min_iter_conv,err{2}(1:min_iter_conv),'-k','DisplayName',num2str(rho_vec(2)));
semilogy(1:min_iter_conv,err{3}(1:min_iter_conv),'-b','DisplayName',num2str(rho_vec(3)));
semilogy(1:min_iter_conv,err{4}(1:min_iter_conv),'-r','DisplayName',num2str(rho_vec(4)));
grid off
leg = legend('AutoUpdate','on','Location','best');
title(leg,'$\rho$');
xlabel('Iteration $j$', 'Interpreter', 'latex') 


rmpath ../
rmpath ../solvers
rmpath ../solvers/expipg
rmpath ../utils