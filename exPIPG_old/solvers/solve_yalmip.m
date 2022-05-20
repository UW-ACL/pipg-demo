function sol = solve_yalmip(pbm,flag)
    yalmip clear
    opts = pbm.yalmip_options;    
    cnstr = [];
    switch flag
        case 'vec'
            Z = sdpvar(pbm.vec_pbm.len_Z,1);
            cnstr = [cnstr;pbm.vec_pbm.Hineq*Z <= pbm.vec_pbm.hineq];
            cnstr = [cnstr;pbm.vec_pbm.Heq*Z == pbm.vec_pbm.heq];
            yal_out = optimize(cnstr,0.5*Z'*pbm.vec_pbm.P*Z+pbm.vec_pbm.p'*Z,pbm.yalmip_options);        
            sol = devec_solution(pbm,value(Z));
            sol.parse_time = yal_out.yalmiptime*1000;
            sol.solve_time = yal_out.solvertime*1000 + sol.parse_time; 
            sol.name = horzcat('YALMIP-',upper(pbm.yalmip_options.solver),' vec.');
            sol.color = [0.3 0.7 0.9];
            strflag = 'vec.  ';
        case 'devec'
            N = pbm.N;
            M = pbm.M;
            A = pbm.Ad;
            B = pbm.Bd;
            dt = pbm.dt;
            qmax = pbm.qmax;
            umax = pbm.umax;
            Q = pbm.Q;
            R = pbm.R;
            qq = pbm.qq; % nx x N
            ru = pbm.ru; % nu x N-1
            q0 = pbm.q0;
            qf = pbm.qf;
        
            q = sdpvar(pbm.nx,N);
            u = sdpvar(pbm.nu,N-1);
            obj = 0;
            for j = 1:N-1
                cnstr = [cnstr;q(:,j+1) == A*q(:,j) + B*u(:,j)];
                cnstr = [cnstr;u(:,j) <= umax];
                cnstr = [cnstr;-umax <= u(:,j)];
                if j>=M
                    cnstr = [cnstr;q(:,j) == qf];
                elseif j > 1
                    cnstr = [cnstr;q(:,j) <= qmax];
                    cnstr = [cnstr;-qmax <= q(:,j)];
                end
                obj = obj + 0.5*q(:,j)'*Q*q(:,j) + qq(:,j)'*q(:,j) + 0.5*u(:,j)'*R*u(:,j) + ru(:,j)'*u(:,j);
            end
            obj = obj + q(:,N)'*Q*q(:,N) + + qq(:,N)'*q(:,N);
            cnstr = [cnstr;q(:,1) == q0];
            cnstr = [cnstr;q(:,N) == qf];    
            % yal_out = optimize(cnstr,0.5*sum(sum( R * (u .^ 2) + 2*ru .* u )) + 0.5*sum(sum( Q * (q .^ 2) + 2*qq .* q )),opts);
            yal_out = optimize(cnstr,obj,opts);
            sol = struct;
            sol.t = 0:dt:(N-1)*dt;
            sol.q = value(q);
            sol.u = value(u);
            sol.Z = vec_solution(pbm,sol.q,sol.u);
            sol.parse_time = yal_out.yalmiptime*1000;
            sol.solve_time = yal_out.solvertime*1000 + sol.parse_time;
            sol.cost = 0.5*sum(sum( R * (value(u) .^ 2) + 2*ru .* value(u) )) + 0.5*sum(sum( Q * (value(q) .^ 2) + 2*qq .* value(q) ));     
            sol.name = horzcat('YALMIP-',upper(pbm.yalmip_options.solver),' devec.');
            sol.color = [0.9 0.7 0.3];
            strflag = 'devec.';
        otherwise
            error('Invalid flag.')
    end
    switch yal_out.problem
        case 0 
            solve_status = '           Solved';
        case 1
            solve_status = 'Primal infeasible';
            sol.cost = inf;
        case 2
            solve_status = ' Primal unbounded';
            sol.cost = -inf;
        case {12,15} % Infeasible or Unbounded
            solve_status = 'Prim infeas/unbnd';
            sol.cost = nan;
        case 3
            solve_status = '   Max iterations';
    end
    % fprintf('YALMIP-%s %s         %s | Run time: %05.1f ms | Cost: %.3f | Solve time: %.1f ms \n',upper(pbm.yalmip_options.solver),strflag,solve_status,sol.parse_time+sol.solve_time,sol.cost,sol.solve_time);
    fprintf('YALMIP-%s %s         %s | Run time: %05.1f ms | Cost: %.3f\n',upper(pbm.yalmip_options.solver),strflag,solve_status,sol.solve_time,sol.cost);
end