function sol = solve_gurobi(pbm)
    model = struct;
    model.Q = pbm.vec_pbm.P/2;
    model.obj = pbm.vec_pbm.p;
    model.A = pbm.vec_pbm.Hdyn;
    model.rhs = pbm.vec_pbm.hdyn;
    model.sense = '=';
    model.lb = pbm.vec_pbm.Zl;
    model.ub = pbm.vec_pbm.Zu;
    params = struct;
    params.OptimalityTol = pbm.gurobi_eps;
    params.FeasibilityTol = pbm.gurobi_eps;
    params.OutputFlag = 0;
    tic
    result = gurobi(model,params);
    solve_time = toc;
    if strcmp(result.status,'OPTIMAL')
        sol = devec_solution(pbm,result.x);
        sol.solve_status = '           Solved';
        % sol.solve_time = result.runtime*1000;
        sol.solve_time = solve_time*1000;        
    elseif strcmp(result.status,'INFEASIBLE') || strcmp(result.status,'INF_OR_UNBD') || strcmp(result.status,'NUMERIC') || strcmp(result.status,'SUBOPTIMAL')
        sol = struct;
        sol.solve_status = 'Primal infeasible';
        sol.solve_time = solve_time*1000;  
        sol.cost = inf;
    elseif strcmp(result.status,'UNBOUNDED')
        sol = struct;
        sol.solve_status = '  Dual infeasible';
        sol.solve_time = solve_time*1000;  
        sol.cost = -inf;
    elseif strcmp(result.status,'ITERATION_LIMIT') 
        sol = devec_solution(pbm,result.x);
        sol.solve_status = '   Max iterations';
        % sol.solve_time = result.runtime*1000;
        sol.solve_time = solve_time*1000;
    else
        error('Unsupported Gurobi status: %s.',result.status);
    end    
    sol.color = [0.6,0.3,0.1];
    sol.name = 'GUROBI';
    fprintf('GUROBI                       %s | Run time: %05.1f ms | Cost: %.3f\n',sol.solve_status,sol.solve_time,sol.cost);    
end