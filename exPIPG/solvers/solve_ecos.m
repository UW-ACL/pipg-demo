function sol = solve_ecos(pp,ppv)
%{
05/15/2022
Purnanand Elango

Solve the QP obtained by vectorizing the optimal control problem via ECOS

ECOS QP format:
     minimize   0.5 xi^T P xi + q^T xi
    subject to  H xi = g
                xi_min <= xi <= xi_max    

Input:
    Structure of problem data (pp)
    Structure of vectorized problem data (ppv)
Output:
    Structure of solution variables and solver status
%}

    sol = struct;
    sol.name = "ECOS";
    sol.status = "Infeasible";

    ecos_opts = ecosoptimset('ABSTOL',pp.ecos_abs_tol,'RELTOL',pp.ecos_rel_tol,'FEASTOL',pp.ecos_abs_tol,'VERBOSE',0);
    tic
    [xi,~,exit_flag,info,~] = ecosqp(ppv.P,ppv.q,[],[],ppv.H,full(ppv.g),ppv.xi_min,ppv.xi_max,ecos_opts);
    sol.solve_time = toc*1000;
    % sol.solve_time = info.time; 
    % sol.obj_val = info.ecosinfo.pcost
    sol.obj_val = 0.5*xi'*ppv.P*xi;
    if exit_flag == 1
        sol.status = "Feasible";
        sol.xi = xi;
        sol.u = reshape(xi(1:pp.m*pp.N),[pp.m,pp.N]);
        sol.x = reshape(xi(pp.m*pp.N+1:end),[pp.n,pp.N+1]);        
    end

end