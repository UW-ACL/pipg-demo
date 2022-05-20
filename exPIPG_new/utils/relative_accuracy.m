function val = relative_accuracy(s1,s2,varargin)
%{
05/01/2022
Purnanand Elango

Compare two feasible solutions
%}
    if strcmp(s1.status,"Feasible") && strcmp(s2.status,"Feasible")
        % val = norm(s1.xi-s2.xi)/norm(s1.xi);
        val = norm(s1.xi-s2.xi,"inf");
        % val = abs(s1.obj_val-s2.obj_val);        
        if nargin == 2
            fprintf("\nComparison of %s (%.2f ms) and %s (%.2f ms)\n",s1.name,s1.solve_time,s2.name,s2.solve_time);
            fprintf("State error  : %.2e\nControl error: %.2e\nCost error   : %.2e\n",norm(s1.x-s2.x),norm(s1.u-s2.u),norm(s1.obj_val-s2.obj_val));
        elseif nargin == 3 && varargin{1} == 0
            fprintf("%s-%s relative accuracy: %.2e\n",s1.name,s2.name,val);
        else
            error("Incorrect number of input arguments");
        end
        
    % figure
    % subplot(1,2,1)
    % plot(s1.x(1,:),s1.x(2,:),'-b');
    % hold on
    % plot(s1.x(1,1),s1.x(2,1),'ob');
    % plot(s2.x(1,:),s2.x(2,:),'--r');
    % plot(s2.x(1,1),s2.x(2,1),'.r','MarkerSize',18);        
    % title("$x$");
    % subplot(1,2,2)
    % plot(s1.u(1,:),s1.u(2,:),'-b');
    % hold on
    % plot(s1.u(1,1),s1.u(2,1),'ob');
    % plot(s2.u(1,:),s2.u(2,:),'--r');
    % plot(s2.u(1,1),s2.u(2,1),'.r','MarkerSize',18);
    % title("$u$");

    else
        warning("One of the solutions is infeasible.")
    end
end