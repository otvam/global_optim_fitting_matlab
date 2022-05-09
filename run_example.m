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
[var_opt, var_fix, var_err, fct_err, format] = get_param_fitting();

% get the solver structure and options
[cache, optimizer] = get_param_solver();

% solve the fitting problem
[param, optim] = get_solver(var_opt, var_fix, var_err, fct_err, format, cache, optimizer);

end

