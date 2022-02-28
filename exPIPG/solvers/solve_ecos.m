function sol = solve_ecos(pbm)
    ecos_opts = ecosoptimset('ABSTOL',pbm.ecos_eps,'RELTOL',pbm.ecos_eps,'FEASTOL',pbm.ecos_eps,'VERBOSE',0);
    tic
    [Z,~,exitflag,info,~] = ecosqp(pbm.vec_pbm.P,pbm.vec_pbm.p,pbm.vec_pbm.Hineq,pbm.vec_pbm.hineq,pbm.vec_pbm.Heq,pbm.vec_pbm.heq,[],[],ecos_opts);
    solve_time = toc;
    % solve_time = info.time;
    sol = devec_solution(pbm,Z);
    sol.name = 'ECOS';
    sol.color = [0.4,0.1,0.8];
    sol.solve_time = solve_time*1000;
    % sol.run_time = run_time*1000;
    if exitflag == 1
        solve_status = '           Solved';
    elseif exitflag == -2
        solve_status = 'Primal infeasible';
    elseif exitflag == -3
        solve_status = '  Dual infeasible';
    elseif exitflag == 0
        solve_status = '   Max iterations';
    else
        error('Unsupported exitflag: %.0f',exitflag);
    end
    sol.solve_status = solve_status;
    fprintf('ECOS                         %s | Run time: %05.1f ms | Cost: %.3f\n',solve_status,sol.solve_time,sol.cost);            
end