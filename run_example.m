function run_example()
% Example for the global fitting code.
%
%    Fit a model with respect to a dataset.
%    
%
%    Thomas Guillod.
%    2021-2022 - BSD License.

close('all')
addpath('solver')
addpath('param')

% get the solver structure and options
[cache, optimizer] = get_param_solver();

% solve the fitting problem
% [var_opt, var_fix, var_err, fct_err, format] = get_param_fitting();
% [param, optim] = get_solver(var_opt, var_fix, var_err, fct_err, format, cache, optimizer);

% solve the optimization problem
[var_opt, var_fix, var_err, fct_err, format] = get_param_optim();
[param, optim] = get_solver(var_opt, var_fix, var_err, fct_err, format, cache, optimizer);

end

