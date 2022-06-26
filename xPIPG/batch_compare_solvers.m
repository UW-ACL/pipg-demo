%{
05/18/2022
Purnanand Elango

Batch script for obtaining average solve times and infeasibility detection 
times (for a fixed problem size) for high and low solver tolerance parameters
%}

% Feasible problems
compare_solvers_feas_high_revised;
compare_solvers_feas_low_revised;

% Infeasible problems
compare_solvers_infeas_high;
compare_solvers_infeas_low;