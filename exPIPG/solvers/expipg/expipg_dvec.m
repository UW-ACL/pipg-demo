function [x1,u1,w1,exit_flag] = expipg_dvec(x1,u1,w1,...                                          % Main solution variables
                                            x2,u2,w2,x3,u3,w3,err_p1,err_d1,err_d2,...            % Auxiliary variables
                                            lam,omg,rho,k_test,k_max,exit_flag,...                % Algorithm parameters
                                            eps_feas,eps_infeas,...                               % Algorithm tolerances
                                            N,Q,R,A,AT,B,BT,...
                                            x_0,x_N,x_min,x_max,u_min,u_max...
                                            ) %#codegen
%{
05/02/2022
Purnanand Elango

Solve the following trajectory optimization problem via exPIPG:
     minimize   x_N^T Q x_N + \sum_{t=0}^{N-1} x_t^T Q x_t + u_t^T R u_t 
    subject to  x_{t+1} = A x_t + B u_t, t = 0,...,N-1,
                u_min <= u_t <= u_max,   t = 0,...,N-1,
                x_min <= x_t <= x_max,   t = 1,...,N-1,
                x_0 = z_0, x_N = z_N.

Infeasibility detection: new test developed by Yue Yu on 04/30/2020
See docs/docs.md for details

Output:
    exit_flag:
       -1 - Maximum iterations reached; no conclusion
        0 - Primal Infeasible
        1 - Primal Dual Feasible (Solved)

%}
    tic
    %% Power iteration for estimating the spectral norm of H^T H
    sig1 = 100.1;
    sig2 = 200.2;
    while abs(sig2-sig1)/sig1 >= 0.005
        sig2 = sig1;

        for t = 1:N
            w1(:,t) = x1(:,t+1) - A*x1(:,t) - B*u1(:,t);
        end

        x1(:,1) = -AT*w1(:,1);
        u1(:,1) = -BT*w1(:,1);
        for t = 2:N
            x1(:,t) = w1(:,t-1) - AT*w1(:, t);
            u1(:,t) = -BT*w1(:,t);
        end
        x1(:,N+1) = w1(:,N);
    
        % Compact syntax for norm of primal vector: sig1 = norm([x1(:);u1(:)]);
        sig1 = 0;
        for t = 1:N
            sig1 = sig1 + dot(x1(:,t),x1(:,t)) + dot(u1(:,t),u1(:,t));
        end
        sig1 = sqrt(sig1 + dot(x1(:,N+1),x1(:,N+1)));

        x1 = x1 ./ sig1;
        u1 = u1 ./ sig1;
    end
    % Buffer the estimated spectral norm
    sig1 = 1.1*sig1;

    %% exPIPG

    % Compute step sizes
    alf = 2/((lam^2 + 4*omg*sig1)^0.5+lam);
    bet = omg*alf;

    % Initialization
    x1 = x2;
    u1 = u2;
    w1 = w2;
    for k = 1:k_max % Run algorithm for at most k_max iterations 
        % Copy previous iterates
        x2 = x1;
        u2 = u1;
        w2 = w1;

        % Projected gradient with proportional-integral feedback of conic constraint affine term
        x1(:,1) = x_0; % Projection on singleton set
        u1(:,1) = max(u_min, min(u_max, u3(:,1) - alf*(R*u3(:,1) - BT*w3(:,1)))); % Projection on box
        for t = 2:N
            x1(:,t) = max(x_min, min(x_max, x3(:,t) - alf*(Q*x3(:,t) + w3(:,t-1) - AT*w3(:,t)))); % Projection on box
            u1(:,t) = max(u_min, min(u_max, u3(:,t) - alf*(R*u3(:,t) - BT*w3(:,t)))); % Projection on box
        end
        x1(:,N+1) = x_N; % Projection on singleton set

        for t = 1:N
            w1(:,t) = w3(:, t) + bet*(2*x1(:,t+1)-x3(:,t+1) - 2*A*x1(:,t)+A*x3(:,t) - 2*B*u1(:,t)+B*u3(:,t));
        end

        % Extrapolation
        x3 = (1-rho) .* x3 + rho .* x1;
        u3 = (1-rho) .* u3 + rho .* u1;
        w3 = (1-rho) .* w3 + rho .* w1;

        if rem(k,k_test) == 0 % Test stopping criteria only every k_test iterations
            err_p1 = max(max(abs(x2-x3)));
            err_p1 = max(err_p1,max(max(abs(u2-u3))));
            err_d1 = max(max(abs(w2-w3)));

            w2 = w3 - w2;
            x2(:,1) = x_0;
            u2(:,1) = (-BT*w2(:,1) >= 0) .* u_min + (-BT*w2(:,1) < 0) .* u_max;
            for t = 2:N
                x2(:,t) = (w2(:,t-1)-AT*w2(:,t) >= 0) .* x_min + (w2(:,t-1)-AT*w2(:,t) < 0) .* x_max;
                u2(:,t) = (-BT*w2(:,t) >= 0) .* u_min + (-BT*w2(:,t) < 0) .* u_max;
            end
            x2(:,N+1) = x_N;
            err_d2 = 0;
            for t = 1:N
                err_d2 = err_d2 + dot(w2(:,t),x2(:,t+1)-A*x2(:,t)-B*u2(:,t));
            end

            % Test for primal dual feasibility
            if err_p1 <= eps_feas && err_d1 <= eps_feas
                fprintf("\nexPIPG DVEC converged in %.0f iterations: PRIMAL DUAL FEASIBLE\n",k);
                exit_flag = 1; 
                break
            % Test for primal infeasibility
            elseif err_d2 >= -eps_infeas
                fprintf("\nexPIPG DVEC converged in %.0f iterations: PRIMAL INFEASIBLE\n",k);                
                exit_flag = 0; % Primal Infeasible
                break
            end
        end
    end
    solve_time = toc*1000;
    if exit_flag == -1
        fprintf("\nexPIPG DVEC did not converge in %.0f iterations: MAX ITERATIONS REACHED\n",k_max);
    end
    fprintf("Solve time   : %5.1f ms\n",solve_time);
    fprintf("Primal error : %9.2e\nDual error   : %9.2e\nInfeas. cert.: %9.2e\n",err_p1,err_d1,err_d2);  
end