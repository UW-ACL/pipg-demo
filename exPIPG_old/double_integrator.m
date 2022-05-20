function pbm = double_integrator(varargin)

    if nargin > 1
        % varargin = {N,'plot',pbm,sol1,sol2,...}
        if strcmp(varargin{2},'plot')
            plot_solution(varargin{3},varargin{4},varargin(5:end))
        else
            error('Incorrect flag.')
        end
        pbm = [];
        return 
    end

    N = varargin{1};

% state:   q = [x1 v1 x2 v2]
% control: u = [u1 u2] 
    pbm = struct;                                                                               % structure with problem data
    pbm.nx = 4;
    pbm.nu = 2; 
    pbm.dt = 0.2;                                                                               % sampling time [s]
    pbm.N = N;                                                                                  % horizon length
    pbm.M = pbm.N-1; % 38;                                                                      % time instant when system reaches final state and stay until end of horizon  
    
    pbm.pscl = 1.0;
    pbm.vscl = 1.0;
    pbm.uscl = 1.0;
    pbm.sclx = diag([pbm.pscl,pbm.vscl,pbm.pscl,pbm.vscl]);
    pbm.sclu = diag([pbm.uscl,pbm.uscl]);    

    % continuous-time system matrices
    Ac = [0,     1,      0,       0;   
          0,     0,      0,       0;
          0,     0,      0,       1;
          0,     0,      0        0];
    Bc = [0, 0;
          1, 0;
          0, 0;
          0, 1];
   
    [pbm.Ad,pbm.Bd] = get_discrete(Ac,Bc,pbm.dt,'zoh');
    pbm.Ad = (pbm.sclx\pbm.Ad)*pbm.sclx;
    pbm.Bd = (pbm.sclx\pbm.Bd)*pbm.sclu;

%     pbm.Q = pbm.sclx*(1*eye(pbm.nx)/pbm.N)*pbm.sclx;                                         % state quadratic weight matrix
%     pbm.R = pbm.sclu*(1*eye(pbm.nu)/pbm.N)*pbm.sclu;                                         % input quadratic weight matrix

%     pbm.Q = (1*eye(pbm.nx)/pbm.N);                                                             % state quadratic weight matrix
%     pbm.R = (1*eye(pbm.nu)/pbm.N);                                                             % input quadratic weight matrix

    pbm.Q = 1*eye(pbm.nx);                                                             % state quadratic weight matrix
    pbm.R = 1*eye(pbm.nu);                                                             % input quadratic weight matrix

    pbm.qq = 0.0*pbm.sclx*ones(pbm.nx,pbm.N);
    pbm.ru = 0.0*pbm.sclu*ones(pbm.nu,pbm.N-1);    

    gam = 4;
    pbm.q0 = pbm.sclx\[0; -1; 0; 0];                                                           % initial state
    pbm.qf = pbm.sclx\[gam; 0; gam; 0];                                                        % final state

    d = 6;    
    pbm.qmax = pbm.sclx\[d; 1; d; 1];
    pbm.umax = pbm.sclu\[2; 2];                                                                % max abs. val of acceleration [m s^-2]

    %%%% PIPG PARAMETERS %%%%
    pbm.pow_eps             = 1e-3;
    pbm.pipg_eps            = 5e-3; % 5e-3
    pbm.infeas_eps          = 5e-2; % 5e-2   
    pbm.pow_max_iter        = 1.5e4;
    pbm.pipg_max_iter       = 5e4;
    pbm.pipg_rho            = 1.9;
    pbm.pipg_omg            = 100;
    pbm.pipg_test_freq      = 15;
end

function [] = plot_solution(pbm,sol,varargin)

varargin = varargin{1};
amp_fac = 4;
figure

subplot(3,1,1)
legend('AutoUpdate','on','FontSize',11,'Location','eastoutside');
hold on
if nargin>2
    for j = 1:length(varargin)
        plot(varargin{j}.q(1,:),varargin{j}.q(3,:),'DisplayName',varargin{j}.name,'Color',varargin{j}.color,'LineStyle','-','LineWidth',amp_fac*(length(varargin)+1-j+1));
    end
end
plot(sol.q(1,:),sol.q(3,:),'-b','DisplayName',sol.name,'LineWidth',amp_fac);
ylabel('Position [m]');
ylim([-pbm.qmax(3),pbm.qmax(3)]);
xlim([-pbm.qmax(1),pbm.qmax(1)]);
daspect([1 1 1]);

subplot(3,1,2)
hold on
if nargin>2
    for j = 1:length(varargin)
        plot(varargin{j}.q(2,:),varargin{j}.q(4,:),'DisplayName',varargin{j}.name,'Color',varargin{j}.color,'LineStyle','-','LineWidth',amp_fac*(length(varargin)+1-j+1));
    end
end
plot(sol.q(2,:),sol.q(4,:),'-b','DisplayName',sol.name,'LineWidth',amp_fac);
ylabel('Velocity [m]');
ylim([-pbm.qmax(4),pbm.qmax(4)]);
xlim([-pbm.qmax(2),pbm.qmax(2)]);
daspect([1 1 1]);

subplot(3,1,3)
hold on
if nargin>2
    for j = 1:length(varargin)
        plot(varargin{j}.u(1,:),varargin{j}.u(2,:),'DisplayName',varargin{j}.name,'Color',varargin{j}.color,'LineStyle','-','LineWidth',amp_fac*(length(varargin)+1-j+1));
    end
end
plot(sol.u(1,:),sol.u(2,:),'-b','DisplayName',sol.name,'LineWidth',amp_fac);
ylabel('Acceleration [m]');
ylim([-pbm.umax(2),pbm.umax(2)]);
xlim([-pbm.umax(1),pbm.umax(1)]);
daspect([1 1 1]);

set(gcf,'Position',[100,100,900,2700]); % [left bottom width height]

end


       