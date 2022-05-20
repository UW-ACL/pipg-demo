%{
04/30/2022
Purnanand Elango

Solve discrete-time optimal control problem using 
exPIPG, OSQP, SCS, QUADPROG and YALMIP-MOSEK
%}

clear variables
close all
clc

addpath solvers
addpath solvers/expipg/
addpath utils

% x_max = 1; u_max = 0.1;     % Random system
x_max = 1; u_max = 0.5;       % Oscillating masses   

pp = problem_data(24,12,20,x_max,u_max,1);
ppv = construct_vectorized(pp);

% codegen_expipg_varsize(pp,ppv,{'VEC','new_infeas'});
% codegen_expipg_varsize(pp,ppv,{'DVEC','new_infeas'});
% codegen_expipg_varsize(pp,ppv,{'VEC','old_infeas'});
% codegen_expipg_varsize(pp,ppv,{'DVEC','old_infeas'});

% codegen_expipg_constsize(pp,ppv,{'VEC','new_infeas'});
% codegen_expipg_constsize(pp,ppv,{'DVEC','new_infeas'});
% codegen_expipg_constsize(pp,ppv,{'VEC','old_infeas'});
% codegen_expipg_constsize(pp,ppv,{'DVEC','old_infeas'});

sol1 = solve_quadprog(pp,ppv);
sol2 = solve_yalmip(pp);
sol3 = solve_osqp(pp,ppv);
sol4 = solve_scs(pp,ppv);
sol5 = solve_expipg(pp,ppv,{'DVEC','mex','new_infeas'});
sol6 = solve_expipg(pp,ppv,{'VEC','mex','new_infeas'});
sol7 = solve_mosek(pp,ppv);
sol8 = solve_ecos(pp,ppv);

relative_accuracy(sol1,sol2);
relative_accuracy(sol1,sol3);
relative_accuracy(sol1,sol4);
relative_accuracy(sol1,sol5);
relative_accuracy(sol1,sol6);
relative_accuracy(sol1,sol7);
relative_accuracy(sol1,sol8);

constraint_violation(sol1,ppv);
constraint_violation(sol2,ppv);
constraint_violation(sol3,ppv);
constraint_violation(sol4,ppv);
constraint_violation(sol5,ppv);
constraint_violation(sol6,ppv);

rmpath solvers
rmpath solvers/expipg/
rmpath utils