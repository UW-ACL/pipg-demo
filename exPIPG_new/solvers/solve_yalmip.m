function sol = solve_yalmip(pp)
%{
05/01/2022
Purnanand Elango

Model the optimal control problem in YALMIP and solve it using the solver
specified in pp.yalmip_solver

Input:
    Structure of problem data (pp)
Output:
    Structure of solution variables and solver status
%}
    sol = struct;
    sol.name = horzcat('YALMIP-',pp.yalmip_solver);
    sol.status = "Infeasible";

    yalmip clear
    tic
    opts = sdpsettings('solver',pp.yalmip_solver,'verbose',0,...
                       'mosek.MSK_DPAR_INTPNT_QO_TOL_PFEAS',pp.mosek_feas_tol,'mosek.MSK_DPAR_INTPNT_QO_TOL_DFEAS',pp.mosek_feas_tol,...
                       'mosek.MSK_DPAR_INTPNT_QO_TOL_INFEAS',pp.mosek_infeas_tol,'mosek.MSK_DPAR_INTPNT_QO_TOL_REL_GAP',pp.mosek_rel_gap,...
                       'mosek.MSK_DPAR_INTPNT_QO_TOL_NEAR_REL',1);
    x = sdpvar(pp.n,pp.N+1);
    u = sdpvar(pp.m,pp.N);
    % Constraints container
    cnstr = [];

    % Objective
    obj = 0;
    for t = 1:pp.N
        cnstr = [cnstr; x(:,t+1) == pp.A*x(:,t) + pp.B*u(:,t);...
                        u(:,t) <= pp.u_max;...
                        u(:,t) >= -pp.u_max];
        if t>=2
            cnstr = [cnstr; x(:,t) <= pp.x_max;...
                            x(:,t) >= -pp.x_max];
        end
        obj = obj + 0.5*x(:,t)'*pp.Q*x(:,t) + 0.5*u(:,t)'*pp.R*u(:,t);
    end
    obj = obj + 0.5*x(:,pp.N+1)'*pp.Q*x(:,pp.N+1);
    cnstr = [cnstr; x(:,pp.N+1) == pp.zN];
    cnstr = [cnstr; x(:,1) == pp.z0];

    % Call parser and solver
    yalmip_out = optimize(cnstr,obj,opts);
    sol.solve_time = toc*1000;

    % fprintf('YALMIP-%s status: %s\n',opts.solver,yalmiperror(yalmip_out.problem));
    if yalmip_out.problem == 0 % Store solution if the solver converges
        sol.x = value(x);
        sol.u = value(u);
        sol.xi = [reshape(sol.u,[pp.m*pp.N,1]);reshape(sol.x,[pp.n*(pp.N+1),1])]; 
        sol.obj_val = value(obj);
        % sol.solve_time = yalmip_out.solvertime;
        % sol.parse_time = yalmip_out.yalmiptime;
        sol.status = "Feasible";
    end
end