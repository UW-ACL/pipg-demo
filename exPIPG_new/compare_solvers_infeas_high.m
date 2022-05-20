%{
05/14/2022
Purnanand Elango

Compare solvers for determining solution to infeasible problem
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
% codegen_expipg_constsize(pp,ppv,{'DVEC','new_infeas'});
% codegen_expipg_constsize(pp,ppv,{'VEC','new_infeas'});

M_avg = 100;

solve_time = struct;
solve_time.osqp = zeros(1,M_avg);
solve_time.scs = zeros(1,M_avg);
solve_time.mosek = zeros(1,M_avg);
solve_time.ecos = zeros(1,M_avg);
solve_time.expipg = zeros(1,M_avg);

% Set seed for random number generation
rng default

% Infeasibility detection tolerance
eps_infeas = 1e-8; tol_type = 'high';
% eps_infeas = 1e-4; tol_type = 'low';

M_pipg = 0;
for j = 1:M_avg

    pp = problem_data(n,m,N,x_max,u_max,1);
    pp.z0(1:pp.n/2) = pp.z0(1:pp.n/2) + 0.7;  % Push initial state towards the boundary of the box
    ppv = construct_vectorized(pp);

    % OSQP    
    pp.osqp_prim_infeas = eps_infeas;    
    sol_osqp = solve_osqp(pp,ppv);
    assert(sol_osqp.status == "Infeasible");

    % SCS
    pp.scs_infeas_tol = eps_infeas;
    sol_scs = solve_scs(pp,ppv);
    assert(sol_scs.status == "Infeasible");

    % MOSEK
    pp.mosek_infeas_tol = eps_infeas;
    sol_mosek = solve_mosek(pp,ppv);
    assert(sol_mosek.status == "Infeasible");

    % exPIPG
    pp.expipg_infeas_tol = eps_infeas;
    sol_expipg = solve_expipg(pp,ppv,{'DVEC','mex','new_infeas'});
    assert(sol_expipg.status == "Infeasible");

    % ECOS    
    pp.ecos_feas_tol = eps_infeas;
    sol_ecos = solve_ecos(pp,ppv);
    assert(sol_ecos.status == "Infeasible");

    % Solve times
    solve_time.osqp(j) = sol_osqp.solve_time;
    solve_time.scs(j) = sol_scs.solve_time;
    solve_time.mosek(j) = sol_mosek.solve_time;
    solve_time.ecos(j) = sol_ecos.solve_time;
    solve_time.expipg(j) = sol_expipg.solve_time;

    if sol_expipg.status == "Infeasible"
        M_pipg = M_pipg + 1;
    else
        solve_time.expipg(j) = 0;
    end

end

solve_time.osqp_avg = mean(solve_time.osqp);
solve_time.scs_avg = mean(solve_time.scs);
solve_time.mosek_avg = mean(solve_time.mosek);
solve_time.ecos_avg = mean(solve_time.ecos);
solve_time.expipg_avg = sum(solve_time.expipg)/M_pipg;

save(horzcat('results/infeas_',tol_type,'_time_n',num2str(n),'_m',num2str(m),'_N',num2str(N)),'solve_time','M_pipg','M_avg');

fprintf('\n\nn = %d, m = %d, N = %d\n\n',n,m,N);
fprintf('Average detection times\n');
fprintf('OSQP        : %5.1f ms\n',solve_time.osqp_avg);
fprintf('SCS         : %5.1f ms\n',solve_time.scs_avg);
fprintf('MOSEK       : %5.1f ms\n',solve_time.mosek_avg);
fprintf('ECOS        : %5.1f ms\n',solve_time.ecos_avg);
fprintf('%s : %5.1f ms\n',sol_expipg.name,solve_time.expipg_avg);

rmpath utils
rmpath solvers
rmpath solvers/expipg