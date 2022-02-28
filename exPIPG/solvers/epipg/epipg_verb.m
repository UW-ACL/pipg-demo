function [q,u,exitflag,v,err,solve_time] = epipg_verb(q,u,v,err,...                                                                                         % main variables (modified)
                                           eps_pow,eps_pipg,eps_infeas,max_iter_pow,max_iter_pipg,test_freq,rho,omg,lam,N,M,Q,R,qq,ru,A,B,AT,BT,q0,qf,qmax,umax ...              % constants
                                           ,verbose ...       
                                           ) %#codegen
    % epipg which can provide per-iteration information (including cost)
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

    conv_test_primal = 100.0;
    % conv_test_primal_prev = coder.nullcopy(0.0);    
    conv_test_dual = 100.0;
    conv_test_dual_prev = coder.nullcopy(0.0);
    primal_infeas_test = coder.nullcopy(0.0);
    % dual_infeas_test = coder.nullcopy(0.0);
    conv_flag = false;
    cost_val = 0.0;
    exitflag = 0;

    %  0 -> Maximum iterations reached; desired solve or infeasibility detection accuracy not met
    %  1 -> Primal and Dual feasible (solved to desired accuracy)
    % -1 -> Primal infeasible
    % -2 -> Dual infeasible (not available)

    if verbose
        fprintf('\n\n');
        fprintf('----------------------------------------------------------------------------------------------------------------\n');
        fprintf('                    ..:: PIPG v0.1 - Proportional-Integral Projected Gradient Method ::..                       \n');
        fprintf('                                          for Trajectory Optimization                                           \n');        
        fprintf('----------------------------------------------------------------------------------------------------------------\n');
        fprintf('|    Iteration   |   Objective   |    Primal Err.   |     Dual Err.     |    Primal Inf.   |     Dual Inf.     |\n');       
        fprintf('----------------------------------------------------------------------------------------------------------------\n');
    end

    % pipg iterations
    k_err = 0;
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
            cost_val = 0.0;
            for j = 1:N-1
                conv_test_primal = max(conv_test_primal,norm(q(:,j)-qp(:,j),'inf'));
                conv_test_primal = max(conv_test_primal,norm(u(:,j)-up(:,j),'inf'));
                conv_test_dual = max(conv_test_dual,norm(v(:,j)-vp(:,j),'inf'));
                cost_val = cost_val + 0.5*q(:,j)'*Q*q(:,j) + qq(:,j)'*q(:,j) + 0.5*u(:,j)'*R*u(:,j) + ru(:,j)'*u(:,j); 
            end
            cost_val = cost_val + 0.5*q(:,N)'*Q*q(:,N) + qq(:,N)'*q(:,N);
            conv_test_primal = max(conv_test_primal,norm(q(:,N)-qp(:,N),'inf'));
            
            k_err = k_err + 1;
            err(k_err,1) = conv_test_primal;
            err(k_err,2) = conv_test_dual;
            primal_infeas_test =  abs(conv_test_dual/conv_test_dual_prev-1);
            % dual_infeas_test = abs(conv_test_primal/conv_test_primal_prev-1);

            if verbose
                fprintf('|      %04.0f      |    %06.1e    |     %06.1e      |      %06.1e      |     %06.1e      |      %06.1e      |\n',k,cost_val,conv_test_primal,conv_test_dual,primal_infeas_test,dual_infeas_test);
            end
    
            if conv_test_primal <= eps_pipg && conv_test_dual <= eps_pipg
                if verbose
                    fprintf("PRIMAL DUAL FEASIBLE\nPIPG converged in %.0f iterations.\n",k);
                end
                conv_flag = true;
                exitflag = 1;
                break
            elseif primal_infeas_test <= eps_infeas
                conv_flag = true;
                exitflag = -1;
                fprintf("PRIMAL INFEASIBLE\nPIPG converged in %.0f iterations.\n",k);
                break
            % elseif dual_infeas_test <= eps_infeas
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