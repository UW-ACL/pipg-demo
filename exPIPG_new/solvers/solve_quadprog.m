function sol = solve_quadprog(pp,ppv)
%{
05/01/2022
Purnanand Elango

Solve the QP obtained by vectorizing the optimal control problem via MATLAB quadprog

Input:
    Structure of problem data (pp)
    Structure of vectorized problem data (ppv)
Output:
    Structure of solution variables and solver status
%}

    sol = struct;
    sol.name = "QUADPROG";
    sol.status = "Infeasible";

    tic
    opts = optimoptions('quadprog','Algorithm','interior-point-convex',...
                        'OptimalityTolerance',pp.quadprog_opt_tol,'ConstraintTolerance',pp.quadprog_feas_tol,...
                        'StepTolerance',pp.quadprog_step_tol,'Display','none');
    [xi,sol.obj_val,exit_flag] = quadprog(ppv.P,[],[],[],ppv.H,ppv.g,ppv.xi_min,ppv.xi_max,ppv.xi{1},opts);
    sol.solve_time = toc*1000;
    
    if exit_flag == 1 % Store solution if solver converges
        sol.xi = xi;
        sol.u = reshape(xi(1:pp.m*pp.N),[pp.m,pp.N]);
        sol.x = reshape(xi(pp.m*pp.N+1:end),[pp.n,pp.N+1]);
        sol.status = "Feasible"; 
    end
end