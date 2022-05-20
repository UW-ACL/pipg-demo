function sol = solve_epipg(pbm,verbose_flag,m_or_mex)

    if strcmp(m_or_mex,'mex')
        name_suffix = '_mex';
    elseif strcmp(m_or_mex,'m')
        name_suffix = '';
    end

    if verbose_flag
        fname = horzcat('epipg_verb',name_suffix);
        [q,u,exitflag,~,~,solve_time] = feval(fname,pbm.qinit,pbm.uinit,pbm.vinit,pbm.err,...                                                                                                                     % main variables (modified)
                      pbm.eps_pow,pbm.eps_pipg,pbm.eps_infeas,pbm.max_iter_pow,pbm.max_iter,pbm.test_freq,pbm.rho,pbm.omg,pbm.lam,pbm.N,pbm.M,pbm.Q,pbm.R,pbm.qq,pbm.ru,pbm.A,pbm.B,pbm.AT,pbm.BT,pbm.q0,pbm.qf,pbm.qmax,pbm.umax ...    % constants
                      ,verbose_flag);
    else
        fname = horzcat('epipg',name_suffix);
        [q,u,exitflag,~,solve_time] = feval(fname,pbm.qinit,pbm.uinit,pbm.vinit,...                                                                                                                     % main variables (modified)
                      pbm.eps_pow,pbm.eps_pipg,pbm.eps_infeas,pbm.max_iter_pow,pbm.max_iter,pbm.test_freq,pbm.rho,pbm.omg,pbm.lam,pbm.N,pbm.M,pbm.Q,pbm.R,pbm.qq,pbm.ru,pbm.A,pbm.B,pbm.AT,pbm.BT,pbm.q0,pbm.qf,pbm.qmax,pbm.umax ...    % constants
                      );    
    end

    sol = struct;
    sol.q = q;
    sol.u = u;
    sol.Z = [reshape([q(:,1:pbm.N-1);u],[(pbm.nx+pbm.nu)*(pbm.N-1),1]);q(:,pbm.N)];
    sol.cost = 0;
    for j = 1:pbm.N-1
        sol.cost = sol.cost + 0.5*q(:,j)'*pbm.Q*q(:,j) + pbm.qq(:,j)'*q(:,j) + 0.5*u(:,j)'*pbm.R*u(:,j) + pbm.ru(:,j)'*u(:,j);
    end
    sol.cost = sol.cost + 0.5*q(:,pbm.N)'*pbm.Q*q(:,pbm.N) + pbm.qq(:,pbm.N)'*q(:,pbm.N);    
    
    switch exitflag
        case 1
            sol.solve_status = '           Solved';
        case 0
            sol.solve_status = '   Max iterations';
        case -1
            sol.solve_status = 'Primal infeasible';
        % case -2
        %    sol.solve_status = '  Dual infeasible';            
    end
    
    sol.name = 'ePIPG';
    sol.solve_time = solve_time*1000;
    sol.color = [0.4,0.3,0.55];
    fprintf('ePIPG                        %s | Run time: %05.1f ms | Cost: %.3f\n',sol.solve_status,sol.solve_time,sol.cost);
end