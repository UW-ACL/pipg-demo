function sol = solve_scs(pbm)
    settings = struct('eps_abs',pbm.scs_eps,'eps_rel',pbm.scs_eps,...
        'acceleration_lookback',0,...
        'verbose',0);
    data = struct;
    data.P = pbm.vec_pbm.P;
    data.A = [pbm.vec_pbm.Heq;pbm.vec_pbm.Hineq];
    data.c = (pbm.vec_pbm.p);
    data.b = ([pbm.vec_pbm.heq;pbm.vec_pbm.hineq]);
    cone.l = size(pbm.vec_pbm.hineq,1);
    cone.z = size(pbm.vec_pbm.heq,1);
    tic
    [x,~,~,info] = scs(data,cone,settings);
    solve_time = toc;
    sol = devec_solution(pbm,x);
    sol.solve_time = solve_time*1000;
    sol.name = 'SCS';
    sol.color = [0.4,0.8,0.9];
    sol.solve_time = info.solve_time;
    exitflag = info.status_val;
    if exitflag == 1
        solve_status = '           Solved';
    elseif exitflag == 2
        solve_status = '       Inaccurate';
    elseif exitflag == -2 || exitflag == -7
        solve_status = 'Primal infeasible';     
    elseif exitflag == -1 || exitflag == -6
        solve_status = ' Primal unbounded';
    else 
        error('Unsupported exitflag: %.0f',exitflag);
    end

    
    sol.solve_status = solve_status;
    fprintf('SCS                          %s | Run time: %05.1f ms | Cost: %.3f\n',solve_status,sol.solve_time,sol.cost);    
end