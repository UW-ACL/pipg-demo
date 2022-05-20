function [pbm,vec_pbm] = problem_data(varargin)
% This function parses the problem data for a quasi-convex minimum-time
% control problem into vectors and matrices that QP solvers accept.

if nargin == 2
    pbm = feval(varargin{1},varargin{2}); % sys_func_name, N
else
    pbm = varargin{1}; % pbm  already constructed
end

nx = pbm.nx;
nu = pbm.nu;
nxnu = nx+nu;
N = pbm.N;
M = pbm.M;

assert(pbm.M<N && pbm.M>1,'Target reaching instant is invalid');

% vectorized quantities
len_Z = nx*N+nu*(N-1);
P = zeros(len_Z,len_Z);
p = zeros(len_Z,1);

% upper bound and lower bound on decision vector for implementing box
% constraint
Zu = [kron(ones(N-1,1),[pbm.qmax;pbm.umax]);pbm.qmax];
Zl = [kron(ones(N-1,1),-[pbm.qmax;pbm.umax]);-pbm.qmax];

% Hfull, hfull - contains all equality and inequality constraints including
% the boundary conditions, dynamics, state and input box constraints
% represented as conic constraint Hfull*Z - hfull \in K where K is carterian product of R^n_+ and {0} 

Hfull                                                                   = zeros(  nx*(N-1) + ...          % dynamics         *equality*   constraint 
                                                                                  2*nxnu*(N-1) + ...      % state input box  *inequality* constraint
                                                                                  nx*(N-M)  + ...         % stay at target   *equality*   constraint
                                                                                  2*nx        ...         % boundary cond.   *equality*   constraint
                                                                                          ,len_Z);
hfull                                                                   = [zeros(nx*(N-1),1);
                                                                           kron(ones(2*N-2,1),[pbm.qmax;pbm.umax]);
                                                                           kron(ones(N-M+1,1),[pbm.qf]);
                                                                           pbm.q0];

for t = 1:N-1
    P((t-1)*nxnu+1:t*nxnu,(t-1)*nxnu+1:t*nxnu)                          = diag([diag(pbm.Q);diag(pbm.R)]);
    p((t-1)*nxnu+1:t*nxnu)                                              = [pbm.qq(:,t);pbm.ru(:,t)];    
    Hfull((t-1)*nx+1:t*nx,(t-1)*nxnu+1:(t-1)*nxnu+nxnu+nx)              = [-pbm.Ad -pbm.Bd eye(nx)];
    Hfull((N-1)*nx+(t-1)*nxnu+1:(N-1)*nx+t*nxnu,(t-1)*nxnu+1:t*nxnu)    = eye(nxnu); 
    Hfull((N-1)*nx+nxnu*(N-1)+(t-1)*nxnu+1:...
      (N-1)*nx+nxnu*(N-1)+t*nxnu,(t-1)*nxnu+1:t*nxnu)                   = -eye(nxnu);
    if t>=M
        Hfull((N-1)*nx+2*nxnu*(N-1)+(t-1-M+1)*nx+1:...
            (N-1)*nx+2*nxnu*(N-1)+(t-M+1)*nx,(t-1)*nxnu+1:t*nxnu)       = [eye(nx) zeros(nx,nu)];
        Zl((t-1)*nxnu+1:(t-1)*nxnu+nx)                                  = pbm.qf;
        Zu((t-1)*nxnu+1:(t-1)*nxnu+nx)                                  = pbm.qf;
    end
end
P((N-1)*nxnu+1:(N-1)*nxnu+nx,(N-1)*nxnu+1:(N-1)*nxnu+nx)                = pbm.Q;
p((N-1)*nxnu+1:(N-1)*nxnu+nx)                                           = pbm.qq(:,N);
Hfull(end-2*nx+1:end-nx,end-nx+1:end)                                   = eye(nx);
Hfull(end-nx+1:end,1:nx)                                                = eye(nx);
Zl((N-1)*nxnu+1:(N-1)*nxnu+nx)                                          = pbm.qf;
Zu((N-1)*nxnu+1:(N-1)*nxnu+nx)                                          = pbm.qf;
Zl(1:nx)                                                                = pbm.q0;
Zu(1:nx)                                                                = pbm.q0;

% only the dynamics constraint Hdyn*Z = hdyn
Hdyn = Hfull(1:nx*(N-1),:);
hdyn = hfull(1:nx*(N-1),1);

% all equality constraints Heq*Z = heq
Heq = Hfull([1:nx*(N-1),end-nx*(N-M+2)+1:end],:);
heq = hfull([1:nx*(N-1),end-nx*(N-M+2)+1:end],1);

% all ineqality constraints formulated as Hineq*Z <= hineq
Hineq = Hfull(nx*(N-1)+1:end-nx*(N-M+2),:);
hineq = hfull(nx*(N-1)+1:end-nx*(N-M+2),1);

% all constraints formulated as hl <= Hlu*Z <= hu
Hlu = Hfull([1:nx*(N-1)+nxnu*(N-1),end-nx*(N-M+2)+1:end],:);
hu = hfull([1:nx*(N-1)+nxnu*(N-1),end-nx*(N-M+2)+1:end],1);
hl = hu; hl(nx*(N-1)+1:nx*(N-1)+nxnu*(N-1)) = -hu(nx*(N-1)+1:nx*(N-1)+nxnu*(N-1));

vec_pbm = struct;
vec_pbm.len_Z = len_Z;
vec_pbm.P = sparse(P);
vec_pbm.p = p;
vec_pbm.Hfull = sparse(Hfull);
vec_pbm.hfull = hfull;
vec_pbm.Hlu = sparse(Hlu);
vec_pbm.hl = hl;
vec_pbm.hu = hu;
vec_pbm.Hdyn = sparse(Hdyn);
vec_pbm.hdyn = hdyn;
vec_pbm.Hineq = sparse(Hineq);
vec_pbm.hineq = hineq;
vec_pbm.Heq = sparse(Heq);
vec_pbm.heq = heq;
vec_pbm.Zl = Zl;
vec_pbm.Zu = Zu;
vec_pbm.Zinit = randn(len_Z,1); vec_pbm.Zinit = vec_pbm.Zinit ./ norm(vec_pbm.Zinit);
pbm.vec_pbm = vec_pbm;

%%% PIPG VEC. PARAMETERS
% pipg_pbm = struct;
% pipg_pbm.eps_pow                    = pbm.pow_eps;
% pipg_pbm.max_iter_pow               = pbm.pow_max_iter;
% pipg_pbm.eps_pipg                   = pbm.pipg_eps;
% pipg_pbm.eps_infeas                 = pbm.infeas_eps;
% pipg_pbm.max_iter                   = pbm.pipg_max_iter; 
% pipg_pbm.test_freq                  = pbm.pipg_test_freq;
% pipg_pbm.rho                        = pbm.pipg_rho;
% pipg_pbm.omg                        = pbm.pipg_omg;
% pipg_pbm.lam                        = full(max(diag(vec_pbm.P)));
% pipg_pbm.P                          = vec_pbm.P;    
% pipg_pbm.p                          = vec_pbm.p;
% pipg_pbm.H                          = vec_pbm.Heq(1:nx*(N-1),:);
% pipg_pbm.HT                         = vec_pbm.Heq(1:nx*(N-1),:)';
% pipg_pbm.h                          = vec_pbm.heq(1:nx*(N-1),1);
% pipg_pbm.Zu                         = vec_pbm.Zu;
% pipg_pbm.Zl                         = vec_pbm.Zl;
% pipg_pbm.Zinit                      = randn(len_Z,1); pipg_pbm.Zinit = pipg_pbm.Zinit ./ norm(pipg_pbm.Zinit);
% pipg_pbm.Winit                      = randn(nx*(N-1),1);
% pipg_pbm.nx                         = nx;
% pipg_pbm.nu                         = nu;
% pipg_pbm.N                          = N;
% pipg_pbm.dt                         = pbm.dt;
% pbm.pipg_pbm = pipg_pbm;
% if N<=400
%     svd_HTH = svd(full(pipg_pbm.HT*pipg_pbm.H));
%     fprintf('PIPG: maximum singular value of H^T*H = %.2f\n',max(svd_HTH));
% end
%%%%

%%% PIPG DEVEC. PARAMETERS
pipg_dvec_pbm = struct;
pipg_dvec_pbm.eps_pow                    = pbm.pow_eps;
pipg_dvec_pbm.max_iter_pow               = pbm.pow_max_iter;
pipg_dvec_pbm.eps_pipg                   = pbm.pipg_eps;
pipg_dvec_pbm.eps_infeas                 = pbm.infeas_eps;
pipg_dvec_pbm.max_iter                   = pbm.pipg_max_iter;
pipg_dvec_pbm.test_freq                  = pbm.pipg_test_freq;
pipg_dvec_pbm.rho                        = pbm.pipg_rho;
pipg_dvec_pbm.omg                        = pbm.pipg_omg;
pipg_dvec_pbm.lam                        = full(max(diag(vec_pbm.P)));
pipg_dvec_pbm.nx                         = pbm.nx;
pipg_dvec_pbm.nu                         = pbm.nu;
pipg_dvec_pbm.N                          = pbm.N;
pipg_dvec_pbm.M                          = pbm.M;
pipg_dvec_pbm.Q                          = pbm.Q;    
pipg_dvec_pbm.R                          = pbm.R;
pipg_dvec_pbm.qq                         = pbm.qq;    
pipg_dvec_pbm.ru                         = pbm.ru;
pipg_dvec_pbm.A                          = pbm.Ad;    
pipg_dvec_pbm.B                          = pbm.Bd;
pipg_dvec_pbm.AT                         = pbm.Ad';    
pipg_dvec_pbm.BT                         = pbm.Bd';
pipg_dvec_pbm.q0                         = pbm.q0;
pipg_dvec_pbm.qf                         = pbm.qf;
pipg_dvec_pbm.qmax                       = pbm.qmax;
pipg_dvec_pbm.umax                       = pbm.umax;
qu                                       = reshape(vec_pbm.Zinit(1:nxnu*(N-1),1),[nxnu,N-1]); 
pipg_dvec_pbm.qinit                      = [qu(1:nx,1:N-1),vec_pbm.Zinit(end-nx+1:end)];
pipg_dvec_pbm.uinit                      = qu(nx+1:nxnu,1:N-1); clear qu;
pipg_dvec_pbm.vinit                      = randn(nx,N-1);
pipg_dvec_pbm.err                        = zeros(pbm.pipg_max_iter,2);   
pipg_dvec_pbm.Z                          = vec_pbm.Zinit;
pbm.pipg_dvec_pbm = pipg_dvec_pbm;
%%%%

% solver settings
solve_tol = 1e-8;

pbm.yalmip_options = sdpsettings('verbose',0,'solver',...
                    'gurobi','gurobi.OptimalityTol',solve_tol,'gurobi.FeasibilityTol',solve_tol);
                ... 'mosek','mosek.MSK_DPAR_INTPNT_QO_TOL_PFEAS',solve_tol,'mosek.MSK_DPAR_INTPNT_QO_TOL_DFEAS',solve_tol,'mosek.MSK_DPAR_INTPNT_QO_TOL_INFEAS',solve_tol,'mosek.MSK_DPAR_INTPNT_QO_TOL_REL_GAP',solve_tol);

pbm.osqp_max_iter           = 1e5;
pbm.osqp_eps                = solve_tol;
pbm.quadprog_eps            = solve_tol;
pbm.ecos_eps                = solve_tol;
pbm.mosek_eps               = solve_tol;
pbm.gurobi_eps              = solve_tol;
pbm.scs_eps                 = solve_tol;  
end