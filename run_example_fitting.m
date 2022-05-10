function run_example_fitting()
% Global fitting example.
%
%    Fit a dummy model with a dummy dataset.
%    Demonstrate the code capabilities:
%        - cascaded solvers
%        - variable scaling
%        - caching
%
%    Thomas Guillod.
%    2021-2022 - BSD License.

close('all')
addpath('solver')
addpath('param')

% get the solver structure and options
[cache, optimizer] = get_param_solver();

% get the fitting problem
[var_opt, var_fix, var_err, fct_err, format] = get_param_fitting();

% call the solver
[param, optim] = get_solver(var_opt, var_fix, var_err, fct_err, format, cache, optimizer);

end

