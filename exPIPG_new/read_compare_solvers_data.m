%{
05/17/2022
Purnanand Elango

Read the solve time and infeasibility detection time statistics for
different solvers
%}

clear variables
close all
clc

n = 256;
m = n/2;
N = 20;

% prob_type = "feas"; verb_type = "solve";
prob_type = "infeas"; verb_type = "detection";

% tol_type = "low";
tol_type = "high";

file_name = "results/" + prob_type + "_" + tol_type + "_time_n" + num2str(n) + "_m" + num2str(m) + "_N" + num2str(N);
load(file_name);

fprintf('\n\nn = %d, m = %d, N = %d\n\n',n,m,N);
fprintf('Average %s times\n',verb_type);
fprintf('OSQP        : %5.1f ms\n',solve_time.osqp_avg);
fprintf('SCS         : %5.1f ms\n',solve_time.scs_avg);
fprintf('MOSEK       : %5.1f ms\n',solve_time.mosek_avg);
fprintf('ECOS        : %5.1f ms\n',solve_time.ecos_avg);
fprintf('exPIPG DVEC : %5.1f ms\n',solve_time.expipg_avg);

if prob_type == "feas"
    fprintf('\n\nAverage constraint violation\n');
    fprintf('OSQP        : %7.1e\n',constr_viol.osqp_avg);
    fprintf('SCS         : %7.1e\n',constr_viol.scs_avg);
    fprintf('MOSEK       : %7.1e\n',constr_viol.mosek_avg);
    fprintf('ECOS        : %7.1e\n',constr_viol.ecos_avg);
    fprintf('exPIPG DVEC : %7.1e\n',constr_viol.expipg_avg);
end