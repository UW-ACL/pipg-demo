function sol = solve_expipg(pp,ppv,alg_type)
%{
05/02/2022
Purnanand Elango   

Input:
    alg_type{1}:
        Vectorized ('VEC')
        De-vectorized ('DVEC')
    alg_type{2}:
        MEX ('mex')
        MATLAB ('m')
    alg_type{3}:
        Old infeasibility detection criteria based on ratio of dual error ('old_infeas')
        New infeasibility detection criteria developed by Yue Yu on 04/30/2022 ('new_infeas')
        New infeasibility and feasibility detection criteria developed by Yue Yu on 04/30/2022 and 06/13/2022 respectively ('new')
%}
    sol = struct;
    sol.name = horzcat('exPIPG ',alg_type{1});
    sol.status = "Infeasible";

    switch alg_type{3}
        case 'old_infeas'
            infeas_tag = '_v2';
        case 'new_infeas'
            infeas_tag = '';
        case 'new'
            infeas_tag = '_v3';    
        otherwise
            error('Invalid infeasibility detection criteria.');
    end

    switch alg_type{2}
        case 'mex'
            func_name = horzcat('expipg_',lower(alg_type{1}),infeas_tag);
            % Given the problem size in pp, construct the name of the file
            % which may exist for solving problems of this size
            func_name_const = horzcat(func_name,'_n',num2str(pp.n),'_m',num2str(pp.m),...
                                      '_N',num2str(pp.N),'_mex');
            % Check if such a file exists
            if exist(func_name_const,"file")
                func_name = func_name_const;
            else % If not, use the file which solves problem of varying size
                func_name = horzcat(func_name,'_mex');
            end
        case 'm'
            func_name = horzcat('expipg_',lower(alg_type{1}),infeas_tag);
    end
    
    switch alg_type{1}
        case 'VEC'
            % Initialize exPIPG inputs
            xi1 = ppv.xi{1};  xi2 = ppv.xi{2};  xi3 = ppv.xi{3};
            eta1 = ppv.eta{1};  w2 = ppv.eta{2};  w3 = ppv.eta{3};
            err_p1 = 100.1; err_d1 = 200.2; err_d2 = 300.3;
            lam = ppv.max_eig_P; omg = pp.expipg_omg; rho = pp.expipg_rho;
            k_test = pp.expipg_test_iter; k_max = pp.expipg_max_iter;
            exit_flag = -1; eps_feas = pp.expipg_feas_tol; eps_infeas = pp.expipg_infeas_tol;
            
            tic
            % Call exPIPG
            [xi1,eta1,exit_flag] = feval(func_name,xi1,eta1,...                                % Main solution variables
                                      xi2,w2,xi3,w3,err_p1,err_d1,err_d2,...                % Auxiliary variables
                                      lam,omg,rho,k_test,k_max,exit_flag,...                % Algorithm parameters
                                      eps_feas,eps_infeas,...                               % Algorithm tolerances
                                      ppv.P,ppv.H,ppv.HT,...                                % Problem constants
                                      ppv.xi_min,ppv.xi_max);
            sol.solve_time = toc*1000;
        case 'DVEC'
            % Initialize exPIPG inputs
            x1 = ppv.x{1};  x2 = ppv.x{2};  x3 = ppv.x{3};
            u1 = ppv.u{1};  u2 = ppv.u{2};  u3 = ppv.u{3};
            w1 = ppv.w{1};  w2 = ppv.w{2};  w3 = ppv.w{3};
            err_p1 = 100.1; err_d1 = 200.2; err_d2 = 300.3;
            lam = ppv.max_eig_P; omg = pp.expipg_omg; rho = pp.expipg_rho;
            k_test = pp.expipg_test_iter; k_max = pp.expipg_max_iter;
            exit_flag = -1; eps_feas = pp.expipg_feas_tol; eps_infeas = pp.expipg_infeas_tol;
            
            tic
            % Call exPIPG
            [x1,u1,w1,exit_flag] = feval(func_name,x1,u1,w1,...                              % Main solution variables
                                        x2,u2,w2,x3,u3,w3,err_p1,err_d1,err_d2,...          % Auxiliary variables
                                        lam,omg,rho,k_test,k_max,exit_flag,...              % Algorithm parameters
                                        eps_feas,eps_infeas,...                             % Algorithm tolerances
                                        pp.N,pp.Q,pp.R,pp.A,pp.A',pp.B,pp.B',...            % Problem constants
                                        pp.z0,pp.zN,-pp.x_max,pp.x_max,-pp.u_max,pp.u_max);

            sol.solve_time = toc*1000;
            xi1 = [reshape(u1,[pp.m*pp.N,1]);reshape(x1,[pp.n*(pp.N+1),1])];    
            eta1 = reshape(-w1,[pp.n*pp.N,1]);
    end
    sol.obj_val = 0.5*xi1'*ppv.P*xi1;
    if exit_flag == 1
        sol.status = "Feasible";
        sol.xi = xi1;
        sol.eta = eta1;
        sol.u = reshape(xi1(1:pp.m*pp.N),[pp.m,pp.N]);
        sol.x = reshape(xi1(pp.m*pp.N+1:end),[pp.n,pp.N+1]);
    end
end