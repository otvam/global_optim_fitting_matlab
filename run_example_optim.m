function run_example_optim()
% Global optimization example.
%
%    Find the minimum of a dummy function.
%    Demonstrate the code capabilities:
%        - cascaded solvers
%        - variable scaling
%        - error metrics
%        - plot / display
%        - caching
%
%    Thomas Guillod.
%    2021-2022 - BSD License.

close('all')
addpath('solver')
addpath('param')

% get the solver structure and options
[cache, optimizer] = get_param_solver();

% get the optimization problem
[var_opt, var_fix, var_err, fct_err, format] = get_param_optim();

% call the solver
[param, optim] = get_solver(var_opt, var_fix, var_err, fct_err, format, cache, optimizer);

end

