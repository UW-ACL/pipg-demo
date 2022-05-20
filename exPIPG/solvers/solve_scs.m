function sol = solve_scs(pp,ppv)
%{
05/01/2022
Purnanand Elango

Solve the QP obtained by vectorizing the optimal control problem via OSQP

SCS Quadratic Cone Problem format:
     minimize   0.5 xi^T P xi + q^T xi
    subject to  A xi + s = b
                s \in K

Box cone is used to represent the box constraints in state and control input
For an example of usage, see: 
https://www.cvxgrp.org/scs/examples/python/mpc.html#py-mpc

Input:
    Structure of problem data (pp)
    Structure of vectorized problem data (ppv)
Output:
    Structure of solution variables and solver status
%}
    sol = struct;
    sol.name = "SCS";
    sol.status  = "Infeasible";

    tic
    scs_settings = struct('eps_abs',pp.scs_abs_tol,'eps_rel',pp.scs_rel_tol,'eps_infeas',pp.scs_infeas_tol,...
                          'acceleration_lookback',0,'verbose',0,'normalize',1,'check_termination',1,...
                          'adaptive_scale',1); % check_termination doesn't work
    scs_data = struct;
    % Vectors c and b should be in dense format
    scs_data.P = ppv.P;
    scs_data.c = zeros((pp.n+pp.m)*pp.N+pp.n,1);
    
    % Select first and last states for imposing boundary constraints
    e_1 = zeros(1,pp.N+1); e_1(1) = 1;
    e_Np1 = zeros(1,pp.N+1); e_Np1(pp.N+1) = 1;    

    % Exclude first and last states while imposing box constraints
    eye_xu = speye((pp.m+pp.n)*pp.N+pp.n);
    eye_xu(pp.m*pp.N+1:pp.m*pp.N+pp.n,:) = [];
    eye_xu(pp.m*pp.N+pp.n*(pp.N-1)+1:pp.m*pp.N+pp.n*(pp.N-1)+pp.n,:) = [];
    
    scs_data.A = [sparse(pp.n,pp.m*pp.N), sparse(kron(e_1,eye(pp.n)));      % Zero cone (initial condition)
                  sparse(pp.n,pp.m*pp.N), sparse(kron(e_Np1,eye(pp.n)));    % Zero cone (final condition)
                  ppv.H;                                                    % Zero cone (dynamics constraints) 
                  sparse(1,(pp.n+pp.m)*pp.N+pp.n);                          % Box cone  (box cone slack variable constrained to 1)
                  -eye_xu];                                                 % Box cone  (box constraint on state and input)
    scs_data.b = [pp.z0;
                  pp.zN;
                  full(ppv.g);
                  1;
                  zeros((pp.n+pp.m)*pp.N-pp.n,1)];

    scs_cone = struct;
    % Dimension of zero cone
    scs_cone.z = pp.n*(pp.N+2); 
    % Box cone upper and lower bound; dimension of cone inferred from length
    scs_cone.bu = [kron(ones(pp.m*pp.N,1),pp.u_max);kron(ones(pp.n*pp.N-pp.n,1),pp.x_max)];  
    scs_cone.bl = [kron(ones(pp.m*pp.N,1),-pp.u_max);kron(ones(pp.n*pp.N-pp.n,1),-pp.x_max)];

    [x,~,~,info] = scs(scs_data,scs_cone,scs_settings);
    exit_flag = info.status_val;
    sol.obj_val = info.pobj;
    sol.solve_time = toc*1000;

    if exit_flag == 1 % Store solution if solver converges
        sol.xi = x;
        sol.u = reshape(x(1:pp.m*pp.N),[pp.m,pp.N]);
        sol.x = reshape(x(pp.m*pp.N+1:end),[pp.n,pp.N+1]);        
        sol.status = "Feasible";
    end
end