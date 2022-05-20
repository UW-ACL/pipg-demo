function sol = solve_mosek(pbm)
    mosek_opts = struct;
    % mosek_opts.Display = 'off'; % doesn't work
    mosek_opts.MSK_DPAR_INTPNT_QO_TOL_PFEAS = pbm.mosek_eps;
    mosek_opts.MSK_DPAR_INTPNT_QO_TOL_DFEAS = pbm.mosek_eps;
    mosek_opts.MSK_DPAR_INTPNT_QO_TOL_INFEAS = pbm.mosek_eps;
    mosek_opts.MSK_DPAR_INTPNT_QO_TOL_REL_GAP = pbm.mosek_eps;
    tic
    res = mskqpopt(pbm.vec_pbm.P,pbm.vec_pbm.p,pbm.vec_pbm.Hlu,pbm.vec_pbm.hl,pbm.vec_pbm.hu,[],[],mosek_opts,'minimize info echo(0)');
    solve_time = toc;
    sol = devec_solution(pbm,res.sol.itr.xx);
    % sol.solve_time = res.info.MSK_DINF_OPTIMIZER_TIME*1000;
    sol.solve_time = solve_time*1000;
    if strcmp(res.sol.itr.solsta,'OPTIMAL')
        sol.solve_status = '           Solved';
    elseif strcmp(res.sol.itr.solsta,'PRIMAL_INFEASIBLE_CER')
        sol.solve_status = 'Primal infeasible';
    elseif strcmp(res.sol.itr.solsta,'DUAL_INFEASIBLE_CER')
        sol.solve_status = '  Dual infeasible'; 
    elseif strcmp(res.sol.itr.solsta,'UNKNOWN')
        sol.solve_status = '  Numerical issue';
        % sol.solve_status = '   Max iterations';
    end
    sol.name = 'MOSEK';
    sol.color = [0 0.8 0.1];
    fprintf('MOSEK                        %s | Run time: %05.1f ms | Cost: %.3f\n',sol.solve_status,sol.solve_time,sol.cost);
end