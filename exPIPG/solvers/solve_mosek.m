function sol = solve_mosek(pp,ppv)
%{
05/14/2022
Purnanand Elango

Solve the QP obtained by vectorizing the optimal control problem via MOSEK

MOSEK QP format:
     minimize   0.5 xi^T P xi + q^T xi
    subject to  g <= H xi <= g
                xi_min <= xi <= xi_max 
Input:
    Structure of problem data (pp)
    Structure of vectorized problem data (ppv)
Output:
    Structure of solution variables and solver status
%}
    sol = struct;
    sol.name = 'MOSEK';
    sol.status = "Infeasible";

    tic    
    % MOSEK setting
    mosek_opts = struct;
    mosek_opts.MSK_DPAR_INTPNT_QO_TOL_PFEAS = pp.mosek_feas_tol;
    mosek_opts.MSK_DPAR_INTPNT_QO_TOL_DFEAS = pp.mosek_feas_tol;
    mosek_opts.MSK_DPAR_INTPNT_QO_TOL_INFEAS = pp.mosek_infeas_tol;
    mosek_opts.MSK_DPAR_INTPNT_QO_TOL_REL_GAP = pp.mosek_rel_gap;
    % Call MOSEK QP solver
    res = mskqpopt(ppv.P,ppv.q,ppv.H,ppv.g,ppv.g,ppv.xi_min,ppv.xi_max,mosek_opts,'minimize info echo(0)');
    sol.obj_val = 0.5*res.sol.itr.xx'*ppv.P*res.sol.itr.xx;
    sol.solve_time = toc*1000;

    if res.sol.itr.solsta == "OPTIMAL"
        sol.xi = res.sol.itr.xx;
        sol.u = reshape(sol.xi(1:pp.m*pp.N),[pp.m,pp.N]);
        sol.x = reshape(sol.xi(pp.m*pp.N+1:end),[pp.n,pp.N+1]);
        sol.status = "Feasible";
    end

end