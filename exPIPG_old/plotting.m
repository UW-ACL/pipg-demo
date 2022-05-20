clear variables
clc
close all

% load data_feasible_100samples
load data_infeasible_100samples

std_scl = 1;
M = length(Nvec);

figure
shadedErrorBar(1:M,stat_gurobi(:,1),std_scl*stat_gurobi(:,2),'lineProps',{'Marker','.','LineStyle','-','Color',[0.7,0,0],'MarkerSize',20,'LineWidth',2,'DisplayName','GUROBI'});
hold on
shadedErrorBar(1:M,stat_mosek(:,1),std_scl*stat_mosek(:,2),'lineProps',{'Marker','.','LineStyle','-','Color',[0,0,0.7],'MarkerSize',20,'LineWidth',2,'DisplayName','MOSEK'});
shadedErrorBar(1:M,stat_ecos(:,1),std_scl*stat_ecos(:,2),'lineProps',{'Marker','.','LineStyle','-','Color',[0.7,0,0.7],'MarkerSize',20,'LineWidth',2,'DisplayName','ECOS'});
shadedErrorBar(1:M,stat_osqp(:,1),std_scl*stat_osqp(:,2),'lineProps',{'Marker','.','LineStyle','-','Color',[0.7,0.7,0],'MarkerSize',20,'LineWidth',2,'DisplayName','OSQP'});
shadedErrorBar(1:M,stat_scs(:,1),std_scl*stat_scs(:,2),'lineProps',{'Marker','.','LineStyle','-','Color',[0,0.7,0.7],'MarkerSize',20,'LineWidth',2,'DisplayName','SCS'});
shadedErrorBar(1:M,stat_epipg(:,1),std_scl*stat_epipg(:,2),'lineProps',{'Marker','.','LineStyle','-','Color',[0.7,0.7,0.7],'MarkerSize',20,'LineWidth',2,'DisplayName','ePIPG'});
ax = gca;
ax.YScale = 'log';
ylabel('[ms]');
xlabel('Horizon length ($N$)');
Nstr = cell(1,M);
for j = 1:M
    Nstr{j} = num2str(Nvec(j));
end
ax.XTickLabel = Nstr;
ax.XTick = 1:M;
ax.XLim = [1-0.1,M+0.1];
ax.YLim = [0.2,2000];
legend('Location','southeast');
if feas_flag
    title({'Feasible problem',horzcat(num2str(Nmc),' samples')});
else
    title({'Infeasible problem',horzcat(num2str(Nmc),' samples')});
end