%{
05/10/2022
Purnanand Elango

Determine termination tolerances of OSQP, SCS, MOSEK and ECOS so that
solutions returned for feasible problems have the same order of magnitude 
constraint violation and relative distance to optimum (from a high accuracy ground truth)

For each solver tolerance parameter value, M_avg problem instances are
solved by the solver to average over the randomness in
    - Initial guess
    - Initial condition

%}

clear variables
close all
clc

addpath ../
addpath ../utils
addpath ../solvers
addpath ../solvers/expipg

% Parameters

log10_eps_low = log10(1e-4);
log10_eps_high= log10(1e-8);

x_max = 1; 
u_max = 0.5;  

n = 4;
m = n/2;
N = 20;

% File name for recording experiment data
case_name = horzcat('../results/match_tol_n',num2str(n),'_','m',num2str(m),'_','N',num2str(N));

M_eps = 20;     % Number of tolerance values
M_avg = 3;     % Number of solves to average over for each tolerance value

% Structure for holding constraint violation data
constr_viol = struct;
constr_viol.osqp = zeros(1,M_eps);
constr_viol.scs = zeros(1,M_eps);
constr_viol.expipg = zeros(1,M_eps);
constr_viol.mosek = zeros(1,M_eps);
constr_viol.ecos = zeros(1,M_eps);

% Structure for holding relative distance to optimum data
rel_opt = struct;
rel_opt.osqp = zeros(1,M_eps);
rel_opt.scs = zeros(1,M_eps);
rel_opt.expipg = zeros(1,M_eps);
rel_opt.mosek = zeros(1,M_eps);
rel_opt.ecos = zeros(1,M_eps);

% Range of logarithm of solver tolerance parameter
eps_range = logspace(-9,0,M_eps);

for j = 1:M_eps
    
    constr_viol.osqp(j) = 0;
    constr_viol.scs(j) = 0;
    constr_viol.expipg(j) = 0;
    constr_viol.mosek(j) = 0;
    constr_viol.ecos(j) = 0;

    rel_opt.osqp(j) = 0;
    rel_opt.scs(j) = 0;
    rel_opt.expipg(j) = 0;
    rel_opt.mosek(j) = 0;
    rel_opt.ecos(j) = 0;  

    M_avg_expipg = 0;
    M_avg_ecos = 0;
    for k = 1:M_avg

        % Generate problem instance
        pp = problem_data(n,m,N,x_max,u_max,1);
        ppv = construct_vectorized(pp);

        % Obtain high-accuracy ground truth
        sol_truth = solve_mosek(pp,ppv);
        assert(sol_truth.status=="Feasible"); % Ensure that the problem instance is feasible
        
        % OSQP
        pp.osqp_abs_tol = eps_range(j);
        pp.osqp_rel_tol = eps;
        pp.osqp_prim_infeas = eps;
        sol_osqp = solve_osqp(pp,ppv);

        if sol_osqp.status == "Feasible"
            constr_viol.osqp(j) = constr_viol.osqp(j) + constraint_violation(sol_osqp,ppv);
            rel_opt.osqp(j) = rel_opt.osqp(j) + relative_accuracy(sol_osqp,sol_truth,0);
        else
            fprintf("%s problem\n",sol_osqp.status);
            error("OSQP could not solve the problem.");
        end
    
        % SCS
        pp.scs_abs_tol = eps_range(j);
        pp.scs_rel_tol = eps;
        pp.scs_infeas_tol = eps;
        sol_scs = solve_scs(pp,ppv);

        if sol_scs.status == "Feasible"
            constr_viol.scs(j) = constr_viol.scs(j) + constraint_violation(sol_scs,ppv);
            rel_opt.scs(j) = rel_opt.scs(j) + relative_accuracy(sol_scs,sol_truth,0);
        else
            fprintf("%s problem\n",sol_scs.status);
            error("SCS could not solve the problem.");
        end       
    
        % ECOS
        pp.ecos_abs_tol = eps_range(j);
        pp.ecos_feas_tol = eps_range(j);
        pp.ecos_rel_tol = eps; 
        sol_ecos = solve_ecos(pp,ppv);

        if sol_ecos.status == "Feasible"
            M_avg_ecos = M_avg_ecos + 1;
            constr_viol.ecos(j) = constr_viol.ecos(j) + constraint_violation(sol_ecos,ppv);
            rel_opt.ecos(j) = rel_opt.ecos(j) + relative_accuracy(sol_ecos,sol_truth,0);
        else
            fprintf("%s problem\n",sol_ecos.status);
            % error("ECOS could not solve the problem.");
        end
    
        % MOSEK
        pp.mosek_feas_tol = eps_range(j);
        pp.mosek_rel_gap = eps_range(j);
        pp.mosek_infeas_tol = eps;    
        sol_mosek = solve_mosek(pp,ppv);

        if sol_mosek.status == "Feasible"
            constr_viol.mosek(j) = constr_viol.mosek(j) + constraint_violation(sol_mosek,ppv);
            rel_opt.mosek(j) = rel_opt.mosek(j) + relative_accuracy(sol_mosek,sol_truth,0);
        else
            fprintf("%s problem\n",sol_mosek.status);
            error("MOSEK could not solve the problem.");
        end
    
        % exPIPG
        pp.expipg_feas_tol = eps_range(j);
        pp.expipg_infeas_tol = eps;
        sol_expipg = solve_expipg(pp,ppv,{'DVEC','mex','new_infeas'});
    
        if sol_expipg.status == "Feasible"
            M_avg_expipg = M_avg_expipg + 1;
            constr_viol.expipg(j) = constr_viol.expipg(j) + constraint_violation(sol_expipg,ppv);
            rel_opt.expipg(j) = rel_opt.expipg(j) + relative_accuracy(sol_expipg,sol_truth,0);
        else
            fprintf("%s problem\n",sol_expipg.status);
        end
   
    end
    constr_viol.osqp(j) = constr_viol.osqp(j)/M_avg;
    constr_viol.scs(j) = constr_viol.scs(j)/M_avg;
    constr_viol.expipg(j) = constr_viol.expipg(j)/M_avg_expipg;
    constr_viol.mosek(j) = constr_viol.mosek(j)/M_avg;

    rel_opt.osqp(j) = rel_opt.osqp(j)/M_avg;
    rel_opt.scs(j) = rel_opt.scs(j)/M_avg;
    rel_opt.expipg(j) = rel_opt.expipg(j)/M_avg_expipg;
    rel_opt.mosek(j) = rel_opt.mosek(j)/M_avg;

    if M_avg_ecos == 0
        constr_viol.ecos(j) = NaN;
        rel_opt.ecos(j) = NaN;
    else
        constr_viol.ecos(j) = constr_viol.ecos(j)/M_avg_ecos;
        rel_opt.ecos(j) = rel_opt.ecos(j)/M_avg_ecos;    
    end

    clear pp ppv

end

figure
subplot(1,2,1)
loglog(eps_range,constr_viol.osqp,'o-b','DisplayName','OSQP');
hold on
loglog(eps_range,constr_viol.scs,'o-r','DisplayName','SCS');
loglog(eps_range,constr_viol.expipg,'o-k','DisplayName','exPIPG');
loglog(eps_range,constr_viol.mosek,'o-m','DisplayName','MOSEK');
loglog(eps_range,constr_viol.ecos,'o-g','DisplayName','ECOS');
legend('AutoUpdate','on','Location','best');
title({'Constraint violation',horzcat('$n =$ ',num2str(n),', $m =$ ',num2str(m),', $N =$ ',num2str(N))});
xlabel('$\epsilon$');

subplot(1,2,2)
loglog(eps_range,rel_opt.osqp,'o-b','DisplayName','OSQP');
hold on
loglog(eps_range,rel_opt.scs,'o-r','DisplayName','SCS');
loglog(eps_range,rel_opt.expipg,'o-k','DisplayName','exPIPG');
loglog(eps_range,rel_opt.mosek,'o-m','DisplayName','MOSEK');
loglog(eps_range,rel_opt.ecos,'o-g','DisplayName','ECOS');
legend('AutoUpdate','on','Location','best');
title({'Relative accuracy',horzcat('$n =$ ',num2str(n),', $m =$ ',num2str(m),', $N =$ ',num2str(N))});
xlabel('$\epsilon$');

savefig(case_name);

% Determine the solver tolerance for desired accuracy
constr_viol.osqp_eps_low = eps_range(get_eps_idx(constr_viol.osqp,log10_eps_low));

constr_viol.scs_eps_low = eps_range(get_eps_idx(constr_viol.scs,log10_eps_low));

constr_viol.mosek_eps_low = eps_range(get_eps_idx(constr_viol.mosek,log10_eps_low));

constr_viol.ecos_eps_low = eps_range(get_eps_idx(constr_viol.ecos,log10_eps_low));

constr_viol.expipg_eps_low = eps_range(get_eps_idx(constr_viol.expipg,log10_eps_low));

constr_viol.osqp_eps_high = eps_range(get_eps_idx(constr_viol.osqp,log10_eps_high));

constr_viol.scs_eps_high = eps_range(get_eps_idx(constr_viol.scs,log10_eps_high));

constr_viol.mosek_eps_high = eps_range(get_eps_idx(constr_viol.mosek,log10_eps_high));

constr_viol.ecos_eps_high = eps_range(get_eps_idx(constr_viol.ecos,log10_eps_high));

constr_viol.expipg_eps_high = eps_range(get_eps_idx(constr_viol.expipg,log10_eps_high));

% Alternate approach which is slightly faulty
% [~,i_low] = min(abs(log10(constr_viol.osqp)-log10_eps_low));
% [~,i_high] = min(abs(log10(constr_viol.osqp)-log10_eps_high));

save(case_name,'constr_viol','rel_opt','eps_range','n','m','N','M_avg','M_eps','log10_eps_high','log10_eps_low');

rmpath ../
rmpath ../utils
rmpath ../solvers
rmpath ../solvers/expipg

% Determine the index of eps_range corresponding to the solver
% tolerance value which yields constraint violation closest to but less
% than exp(log10_eps)
function idx = get_eps_idx(constr_viol,log10_eps)
    M_eps = length(constr_viol);
    for j = M_eps:-1:1
        if log10(constr_viol(j)) < log10_eps
            idx = j;
            break
        end
    end
end