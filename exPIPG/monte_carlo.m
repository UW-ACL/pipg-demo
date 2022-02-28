clear variables
close all
clc

addpath('utils');
addpath('solvers');
addpath('solvers/epipg');
sys_model = 'double_integrator';

%%% feasible problem
feas_flag = true;
% set manually
% ------------
% gam = 4;
% d = 6;
% omg = 100;
% pscl = 1;
% Q = Inx; R = Inu;
% test_freq = 15;
% eps_pipg = 5e-3;

%%% infeasible problem
% feas_flag = false;
% set manually
% ------------
% gam = 8;
% d = 6;
% omg = 100;
% pscl = 1;
% Q = Inx; R = Inu;
% test_freq = 15;
% eps_infeas = 5e-2;

Nvec = [50,100,200,400,800,1600,3200];
M = length(Nvec);
Nmc = 2;

data_gurobi(M,Nmc) = 0;
data_mosek(M,Nmc) = 0;
data_ecos(M,Nmc) = 0;
data_osqp(M,Nmc) = 0;
data_scs(M,Nmc) = 0;
data_epipg(M,Nmc) = 0;


for j = 1:M
    pbm = feval(sys_model,Nvec(j));
    pbm = problem_data(pbm);
    pipg_pbm = pbm.pipg_dvec_pbm;    
    for k = 0:Nmc
        if k == 0
            codegen_epipg(pipg_pbm,0);
        end

        % reset initialization
        Zinit = randn(pbm.vec_pbm.len_Z,1); 
        Zinit = Zinit ./ norm(Zinit); % normalization
        qu = reshape(Zinit(1:(pbm.nx+pbm.nu)*(pbm.N-1),1),[pbm.nx+pbm.nu,pbm.N-1]); 
        pipg_pbm.qinit = [qu(1:pbm.nx,1:pbm.N-1),Zinit(end-pbm.nx+1:end)];
        pipg_pbm.uinit = qu(pbm.nx+1:(pbm.nx+pbm.nu),1:pbm.N-1); clear qu Zinit;
        pipg_pbm.vinit = randn(pbm.nx,pbm.N-1);

        sol_gurobi = solve_gurobi(pbm);
        sol_mosek = solve_mosek(pbm);
        sol_osqp = solve_osqp(pbm);
        sol_ecos = solve_ecos(pbm);
        sol_scs = solve_scs(pbm);
        sol_epipg = solve_epipg(pipg_pbm,false,'mex');
    
        if feas_flag
            print_accuracy(sol_epipg,sol_gurobi); fprintf("\n");        
        end

        if k>0
            data_gurobi(j,k) = sol_gurobi.solve_time;
            data_mosek(j,k) = sol_mosek.solve_time;
            data_ecos(j,k) = sol_ecos.solve_time;
            data_osqp(j,k) = sol_osqp.solve_time;
            data_scs(j,k) = sol_scs.solve_time;
            data_epipg(j,k) = sol_epipg.solve_time;
        end

    end
end
                    % mean                  % stabdard deviation (div by Nmc-1)
stat_gurobi     = [mean(data_gurobi,2), std(data_gurobi,0,2)];
stat_mosek      = [mean(data_mosek,2),  std(data_mosek,0,2)];
stat_ecos       = [mean(data_ecos,2),   std(data_ecos,0,2)];
stat_osqp       = [mean(data_osqp,2),   std(data_osqp,0,2)];
stat_scs        = [mean(data_scs,2),    std(data_scs,0,2)];
stat_epipg      = [mean(data_epipg,2),  std(data_epipg,0,2)];

if feas_flag
    file_name = horzcat('data_feasible_',num2str(Nmc),'samples');
else
    file_name = horzcat('data_infeasible_',num2str(Nmc),'samples');
end

% save(file_name,...
%     'data_epipg','data_scs','data_ecos','data_osqp','data_mosek','data_gurobi',...
%     'stat_epipg','stat_scs','stat_ecos','stat_osqp','stat_mosek','stat_gurobi','Nvec','Nmc','feas_flag');