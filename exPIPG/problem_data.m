function pp = problem_data(n,m,N,x_max,u_max,varargin)
%{
04/30/2022
Purnanand Elango

Set problem parameters
Random system regulation problem

Input:
    varargin{1}:
        Choice of plant (flag)
            0 - Random neutrally stable system
            1 - System of oscillating masses
    varargin{2}:
        Seed for random number generator
%}

    % Set seed in random number generator for repeatability 
    if nargin == 7
        rng(varargin{2});
    elseif nargin > 7
        error("Invalid number of input arguments.");
    end
    
    % Structure of containing problem parameters
    pp = struct;

    % State dimension
    pp.n = n;

    % Control dimension
    pp.m = m;

    % Horizon length
    pp.N = N;

    % Sampling time (only applicable for oscillating masses)
    pp.T = 0.1;
    
    % Bound on state inf norm
    pp.x_max = x_max;

    % Bound on control input inf norm
    pp.u_max = u_max;

    % Random initial state (which satisfies bound on state)
    pp.z0 = max(-pp.x_max, min(pp.x_max, ...
                [0.1*ones(pp.n/2,1);zeros(pp.n/2,1)] + 0.05*randn(pp.n,1) ...
                ));

    % Final state at origin
    pp.zN = zeros(pp.n,1);           

    % Plant definition
    % flag
    %   0 - Random neutrally stable system
    %   1 - System of oscillating masses
    if nargin == 6
        flag = varargin{1};
    else
        flag = 1;
    end
    pp = plant_defn(pp,flag);
    
    % Weighting matrices in the quadratic cost
    pp.Q = eye(pp.n);
    pp.R = eye(pp.m);
        
    % Solver and parser settings
    pp.expipg_feas_tol = 1e-8;      % Absolute tolerance for successive difference in primal dual variables 
    pp.expipg_infeas_tol = 1e-8;    % Relative tolerance for primal infeasibility detection
    pp.expipg_omg = 200;            % Step-size ratio
    pp.expipg_rho = 1.6;            % Extrapolation factor
    pp.expipg_test_iter = 5e1;      % Stopping criteria testing frequency
    pp.expipg_max_iter = 5e4;       % Maximum iterations

    % Reference for parameters 
    pp.yalmip_solver = 'MOSEK';
    
    % https://docs.mosek.com/latest/rmosek/parameters.html#mosek.dparam.intpnt_qo_tol_dfeas
    pp.mosek_feas_tol = 1e-9;
    pp.mosek_infeas_tol = 1e-9;
    pp.mosek_rel_gap = 1e-9;
    
    pp.quadprog_feas_tol = 1e-12;
    pp.quadprog_opt_tol = 1e-12;
    pp.quadprog_step_tol = 1e-12;

    pp.ecos_abs_tol = 1e-9;
    pp.ecos_rel_tol = 1e-9;
    pp.ecos_feas_tol = 1e-9;
    
    pp.osqp_abs_tol = 1e-9;
    pp.osqp_rel_tol = 1e-9;
    pp.osqp_prim_infeas = 1e-9;
    
    pp.scs_abs_tol = 1e-9;
    pp.scs_rel_tol = 1e-9;  
    pp.scs_infeas_tol = 1e-9;

    % Code generation parameters
    pp.n_max = 256;                 % Max state dim
    pp.m_max = 128;                 % Max control dim
    pp.N_max = 30;                  % Max horizon length

    switch flag
        case 0
            fprintf("..:: Random System ::..\n");
            fprintf("-----------------------\n\n");
        case 1
            fprintf("..:: Oscillating Masses ::..\n");
            fprintf("----------------------------\n\n");
    end
end





