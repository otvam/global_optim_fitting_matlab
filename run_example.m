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

% get the parameter definition and the error function
[var_opt, var_fix, fct_err, format] = get_param_fitting();
% [var_opt, var_fix, fct_err, format] = get_param_optim();

% get the solver structure and options
solver = get_param_solver();

% solve the fitting problem
[param, optim] = get_solver(var_opt, var_fix, fct_err, format, solver);

end

