function sol = solve_osqp(pbm)
    osqp_obj = osqp;
    tic
    osqp_obj.setup(pbm.vec_pbm.P,pbm.vec_pbm.p,pbm.vec_pbm.Hlu,pbm.vec_pbm.hl,pbm.vec_pbm.hu,'eps_abs',pbm.osqp_eps,'eps_rel',pbm.osqp_eps,...
        'max_iter',pbm.osqp_max_iter,...
        'verbose',false);
    sol_osqp = osqp_obj.solve();
    solve_time = toc;
    sol = devec_solution(pbm,sol_osqp.x);
    sol.name = 'OSQP';
    sol.color = [0.5 0.6 0.2];
    % sol.solve_time = sol_osqp.info.run_time*1000;
    sol.solve_time = solve_time*1000;
    solve_status = sol_osqp.info.status;
    sol.solve_status = solve_status;
    exitflag = sol_osqp.info.status_val;
    if exitflag == 1
        solve_status = '           Solved';
    elseif exitflag == 2
        solve_status = '       Inaccurate';
    elseif exitflag == -2
        solve_status = '   Max iterations';
    elseif exitflag == -3
        solve_status = 'Primal infeasible';
    elseif exitflag == -4
        solve_status = '  Dual infeasible';    
    end    
    fprintf('OSQP                         %s | Run time: %05.1f ms | Cost: %.3f\n',solve_status,sol.solve_time,sol.cost);        
end