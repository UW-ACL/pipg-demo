function sol = solve_quadprog(pbm)
    quadprog_opts = optimoptions('quadprog','Display','none','ConstraintTolerance',pbm.quadprog_eps,'OptimalityTolerance',pbm.quadprog_eps);
    tic;
    [Z,~,exitflag] = quadprog(pbm.vec_pbm.P,pbm.vec_pbm.p,pbm.vec_pbm.Hineq,pbm.vec_pbm.hineq,pbm.vec_pbm.Heq,pbm.vec_pbm.heq,[],[],[],quadprog_opts);
    solve_time = toc;
    sol = devec_solution(pbm,Z);
    sol.name = 'QUADPROG';
    sol.color = [0.8,0.4,0.8];
    sol.solve_time = solve_time*1000;
    if exitflag == 1
        solve_status = '           Solved';
    elseif exitflag == -2
        solve_status = 'Primal infeasible';
    elseif exitflag == -3
        solve_status = ' Primal unbounded';        
    elseif exitflag == 0
        solve_status = '   Max iterations';
    end
    sol.solve_status = solve_status;
    fprintf('QUADPROG                     %s | Run time: %05.1f ms | Cost: %.3f\n',solve_status,sol.solve_time,sol.cost);
end