function [q,u,exitflag,v,solve_time] = epipg(q,u,v,...                                                                                                                % main variables (modified)
                                           eps_pow,eps_pipg,eps_infeas,max_iter_pow,max_iter_pipg,test_freq,rho,omg,lam,N,M,Q,R,qq,ru,A,B,AT,BT,q0,qf,qmax,umax ...              % constants      
                                           ) %#codegen
    % barebones epipg which is not verbose and doesn't compute cost at each iteration
    % test for dual infeasibility is turned off

    % additional copies of primal and dual variables
    qp = coder.nullcopy(q);
    up = coder.nullcopy(u);
    vp = coder.nullcopy(v);
    qp1 = coder.nullcopy(q);
    up1 = coder.nullcopy(u);
    vp1 = coder.nullcopy(v);
    q2 = coder.nullcopy(q);
    u2 = coder.nullcopy(u);

    % start epipg timing
    tic
    sig = 1.0;
    sig1 = coder.nullcopy(sig);

    % power iteration
    coder.loop.parallelize('never');
    for k = 1:max_iter_pow
        sig1 = sig;
        for t = 1:N-1
            v(:, t) = q(:, t+1) - A*q(:, t) - B*u(:, t);
        end
        
        q(:, 1) = -AT*v(:, 1);
        u(:, 1) = -BT*v(:, 1);
        for t = 2:N-1
            q(:, t) = v(:, t-1) - AT*v(:, t);
            u(:, t) = -BT*v(:, t);
        end
        q(:, N) = v(:, N-1); 
 
        sig = norm([q(:);u(:)]);

        if abs(sig1-sig) < eps_pow
            break
        end
        q = q ./ sig;
        u = u ./ sig;
        
    end
    sig = 1.1*sig;
    
    % pipg step size

    alf = 2/(sqrt(lam^2 + 4*omg*sig)+lam);
    bet = omg*alf;
    
    % pipg diagnostic variables

    conv_test_primal = coder.nullcopy(0.0);
    % conv_test_primal_prev = coder.nullcopy(0.0);    
    conv_test_dual = coder.nullcopy(0.0);
    conv_test_dual_prev = coder.nullcopy(0.0);
    conv_flag = false;
    exitflag = 0;

    %  0 -> Maximum iterations reached; desired solve or infeasibility detection accuracy not met
    %  1 -> Primal and Dual feasible (solved to desired accuracy)
    % -1 -> Primal infeasible
    % -2 -> Dual infeasible (not available)
    
    % pipg iterations
    coder.loop.parallelize('never');
    for k = 1:max_iter_pipg

       % copy of the previous iterate, used for monitoring convergence
        qp1 = qp;
        up1 = up;
        vp1 = vp;
        
        % proj grad: t = 1
        qp(:, 1) = q0; % projection onto initial value set
        up(:, 1) = max(-umax, min( umax, u(:, 1) - alf*(R*u(:, 1) + ru(:, 1) - BT*v(:, 1)))); % proj grad on u
        
        % proj grad: t in [2, M-1]
        for t = 2:M-1
            qp(:, t) = max(-qmax, min(qmax, q(:, t) - alf*(Q*q(:, t) + qq(:, t) + v(:, t-1) - AT*v(:, t))));
            up(:, t) = max(-umax, min(umax, u(:, t) - alf*(R*u(:, t) + ru(:, t) - BT*v(:, t))));
        end
        
        % proj grad t in [M, N-1]
        for t = M:N-1
            qp(:, t) = qf;
            up(:, t) = max(-umax, min(umax, u(:, t) - alf*(R*u(:, t) + ru(:, t) - BT*v(:, t))));
        end
        
        qp(:, N) = qf;
        
        q2 = 2*qp - q;
        u2 = 2*up - u;
        
        for t = 1:N-1
            vp(:, t) = v(:, t) + bet*(q2(:, t+1) - A*q2(:, t) - B*u2(:, t));
        end
        
        q = (1-rho)*q + rho*qp;
        u = (1-rho)*u + rho*up;
        v = (1-rho)*v + rho*vp;
        
        if rem(k,test_freq)==0
            conv_test_dual_prev = conv_test_dual;
            % conv_test_primal_prev = conv_test_primal;
            conv_test_primal = 0.0; 
            conv_test_dual = 0.0;
            for j = 1:N-1
                conv_test_primal = max(conv_test_primal,norm(q(:,j)-qp(:,j),'inf'));
                conv_test_primal = max(conv_test_primal,norm(u(:,j)-up(:,j),'inf'));
                conv_test_dual = max(conv_test_dual,norm(v(:,j)-vp(:,j),'inf'));
            end
            conv_test_primal = max(conv_test_primal,norm(q(:,N)-qp(:,N),'inf'));
    
            if conv_test_primal <= eps_pipg && conv_test_dual <= eps_pipg
                conv_flag = true;
                exitflag = 1;
                fprintf("PRIMAL DUAL FEASIBLE\nPIPG converged in %.0f iterations.\n",k);
                break
            elseif abs(conv_test_dual/conv_test_dual_prev-1) <= eps_infeas
                conv_flag = true;
                exitflag = -1;
                fprintf("PRIMAL INFEASIBLE\nPIPG converged in %.0f iterations.\n",k);
                break
            % elseif abs(conv_test_primal/conv_test_primal_prev-1) <= eps_infeas
            %     conv_flag = true;
            %     exitflag = -2;
            %     fprintf("DUAL INFEASIBLE\nPIPG converged in %.0f iterations.\n",k);
            %     break
            end
        end

    end
    solve_time = toc; % end epipg timing
    if ~conv_flag
        fprintf("INCONCLUSIVE\nMaximum iterations reached.\n");
    end
end