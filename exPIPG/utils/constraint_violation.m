function val = constraint_violation(sol,ppv)
%{
05/10/2022
Purnanand Elango

Return maximum violation in constraints
    H xi - g = 0
    xi_min <= xi <= xi_max

%}

val = max(norm(ppv.H*sol.xi - ppv.g,"inf"),...
          norm(sol.xi - min(ppv.xi_max, max(ppv.xi_min, sol.xi)),"inf"));

fprintf("%s maximum constraint violation: %.2e\n",sol.name,val);

end