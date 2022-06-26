%{
05/14/2022
Purnanand Elango

Compare solvers for determining solution to feasible problem
with coarse termination tolerance
Low accuracy solution is estimated
%}

clc
clear variables
close all

addpath utils
addpath solvers
addpath solvers/expipg

x_max = 1; 
u_max = 0.5;

n = 256;
m = n/2;
N = 20;

% Code generation (not required for certain fixed problem sizes)
% pp = problem_data(n,m,N,x_max,u_max,1);
% ppv = construct_vectorized(pp);
% codegen_expipg_constsize(pp,ppv,{'DVEC','new'});
% codegen_expipg_constsize(pp,ppv,{'VEC','new'});

cv_struct = struct;

M_avg = 100;
eps_high = 1e-8;
log10_eps_high = log10(eps_high);

solve_time = struct;
solve_time.osqp = zeros(1,M_avg);
solve_time.scs = zeros(1,M_avg);
solve_time.mosek = zeros(1,M_avg);
solve_time.ecos = zeros(1,M_avg);
solve_time.expipg = zeros(1,M_avg);

constr_viol = struct;
constr_viol.osqp = zeros(1,M_avg);
constr_viol.scs = zeros(1,M_avg);
constr_viol.mosek = zeros(1,M_avg);
constr_viol.ecos = zeros(1,M_avg);
constr_viol.expipg = zeros(1,M_avg);

M_ecos = 0;
M_pipg = 0;
for j = 1:M_avg

    pp = problem_data(n,m,N,x_max,u_max,1);
    ppv = construct_vectorized(pp);

    % OSQP    
    pp.osqp_abs_tol = eps_high;
    pp.osqp_rel_tol = eps;
    pp.osqp_prim_infeas = eps;    
    sol_osqp = solve_osqp(pp,ppv);
    constr_viol.osqp(j) = constraint_violation(sol_osqp,ppv);

    % SCS
    pp.scs_abs_tol = eps_high;
    pp.scs_rel_tol = eps;
    pp.scs_infeas_tol = eps;
    sol_scs = solve_scs(pp,ppv);
    constr_viol.scs(j) = constraint_violation(sol_scs,ppv);

    % MOSEK
    pp.mosek_feas_tol = eps_high;
    pp.mosek_rel_gap = eps_high;
    pp.mosek_infeas_tol = eps;
    sol_mosek = solve_mosek(pp,ppv);
    constr_viol.mosek(j) = constraint_violation(sol_mosek,ppv);

    % exPIPG
    pp.expipg_feas_tol = eps_high;
    pp.expipg_infeas_tol = eps;
    sol_expipg = solve_expipg(pp,ppv,{'DVEC','mex','new'});

    % ECOS    
    pp.ecos_abs_tol = eps_high;
    pp.ecos_rel_tol = eps;
    pp.ecos_feas_tol = eps_high;
    sol_ecos = solve_ecos(pp,ppv);

    % Solve times
    solve_time.osqp(j) = sol_osqp.solve_time;
    solve_time.scs(j) = sol_scs.solve_time;
    solve_time.mosek(j) = sol_mosek.solve_time;
    solve_time.ecos(j) = sol_ecos.solve_time;
    solve_time.expipg(j) = sol_expipg.solve_time;

    if sol_ecos.status == "Feasible"
        constr_viol.ecos(j) = constraint_violation(sol_ecos,ppv);
        M_ecos = M_ecos + 1;
    else
        constr_viol.ecos(j) = 0;  
        solve_time.ecos(j) = 0;
    end

    if sol_expipg.status == "Feasible"
        constr_viol.expipg(j) = constraint_violation(sol_expipg,ppv);
        M_pipg = M_pipg + 1;
    else
        constr_viol.expipg(j) = 0;
        solve_time.expipg(j) = 0;
    end

end

solve_time.osqp_avg = mean(solve_time.osqp);
solve_time.scs_avg = mean(solve_time.scs);
solve_time.mosek_avg = mean(solve_time.mosek);
solve_time.ecos_avg = sum(solve_time.ecos)/M_ecos;
solve_time.expipg_avg = sum(solve_time.expipg)/M_pipg;

constr_viol.osqp_avg = mean(constr_viol.osqp);
constr_viol.scs_avg = mean(constr_viol.scs);
constr_viol.mosek_avg = mean(constr_viol.mosek);
constr_viol.ecos_avg = sum(constr_viol.ecos)/M_ecos;
constr_viol.expipg_avg = sum(constr_viol.expipg)/M_pipg;

fprintf('\n\nn = %d, m = %d, N = %d\n\n',n,m,N);
fprintf('Average solve times\n');
fprintf('OSQP        : %5.1f ms\n',solve_time.osqp_avg);
fprintf('SCS         : %5.1f ms\n',solve_time.scs_avg);
fprintf('MOSEK       : %5.1f ms\n',solve_time.mosek_avg);
fprintf('ECOS        : %5.1f ms\n',solve_time.ecos_avg);
fprintf('%s : %5.1f ms\n',sol_expipg.name,solve_time.expipg_avg);

fprintf('\n\nAverage constraint violation\n');
fprintf('OSQP        : %7.1e\n',constr_viol.osqp_avg);
fprintf('SCS         : %7.1e\n',constr_viol.scs_avg);
fprintf('MOSEK       : %7.1e\n',constr_viol.mosek_avg);
fprintf('ECOS        : %7.1e\n',constr_viol.ecos_avg);
fprintf('%s : %7.1e\n',sol_expipg.name,constr_viol.expipg_avg);

save(horzcat('results/feas_high_time_n',num2str(n),'_m',num2str(m),'_N',num2str(N),'_revised'),'solve_time','M_pipg','M_avg','log10_eps_high','constr_viol');

rmpath utils
rmpath solvers
rmpath solvers/expipg