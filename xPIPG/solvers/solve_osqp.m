function sol = solve_osqp(pp,ppv)
%{
05/01/2022
Purnanand Elango

Solve the QP obtained by vectorizing the optimal control problem via OSQP

OSQP QP format:
     minimize   0.5 xi^T P xi + q^T xi
    subject to  l <= A xi <= u

Input:
    Structure of problem data (pp)
    Structure of vectorized problem data (ppv)
Output:
    Structure of solution variables and solver status
%}

    sol = struct;
    sol.name = "OSQP";
    sol.status = "Infeasible";

    tic
    % Contruct parameters according to OSQP QP format
    P = ppv.P;
    q = sparse(pp.m*pp.N+pp.n*pp.N+pp.n,1);
    A = [ppv.H;eye(pp.m*pp.N+pp.n*pp.N+pp.n)];
    l = [ppv.g;ppv.xi_min];
    u = [ppv.g;ppv.xi_max];

    % Setup OSQP problem object and solve
    osqp_obj = osqp; 
    osqp_obj.setup(P,q,A,l,u,'eps_abs',pp.osqp_abs_tol,'eps_rel',pp.osqp_rel_tol,...
                   'eps_prim_inf',pp.osqp_prim_infeas,...
                   'verbose',false,...
                   'polish',false,'check_termination',1);
    osqp_sol = osqp_obj.solve();
    exit_flag = osqp_sol.info.status_val;
    sol.obj_val = osqp_sol.info.obj_val;
    sol.solve_time = toc*1000;

    if exit_flag == 1 % Store solution if solver converges
        sol.xi = osqp_sol.x;
        sol.u = reshape(osqp_sol.x(1:pp.m*pp.N),[pp.m,pp.N]);
        sol.x = reshape(osqp_sol.x(pp.m*pp.N+1:end),[pp.n,pp.N+1]);        
        sol.status = "Feasible";
    end
end