function [var_opt, var_fix, var_err, fct_err, format] = get_param_optim()

var_fix = {};
var_fix{end+1} = struct('name', 'log_opt', 'x0', 15.0, 'idx', 1);
var_fix{end+1} = struct('name', 'lin_opt', 'x0', 1.5, 'idx', 1);
var_fix{end+1} = struct('name', 'lin_opt', 'x0', 2.5, 'idx', 2);
var_fix{end+1} = struct('name', 'lin_opt', 'x0', 0.5, 'idx', 3);

var_opt = {};
var_opt{end+1} = struct('name', 'log_var', 'x0', 1.0, 'lb', 0.01, 'ub', 100.0, 'tol_bnd', 0.05, 'trf', 'log', 'norm', true, 'idx', 1);
var_opt{end+1} = struct('name', 'lin_var', 'x0', 2.0, 'lb', 0.0, 'ub', 3.0, 'tol_bnd', 0.05, 'trf', 'lin', 'norm', true, 'idx', 1);
var_opt{end+1} = struct('name', 'lin_var', 'x0', 2.0, 'lb', 0.0, 'ub', 3.0, 'tol_bnd', 0.05, 'trf', 'lin', 'norm', true, 'idx', 2);
var_opt{end+1} = struct('name', 'lin_var', 'x0', 2.0, 'lb', 0.0, 'ub', 3.0, 'tol_bnd', 0.05, 'trf', 'lin', 'norm', true, 'idx', 3);

var_err = struct('type', 'error', 'arg', 'avg');

fct_err = @(param, n) get_fct_err(param, n);

format.err = struct('spec', '%.3f', 'scale', 1e0, 'unit', 'a.u.');
format.param.log_opt = struct('spec', '%.3f', 'scale', 1e0, 'unit', 'a.u.');
format.param.lin_opt = struct('spec', '%.3f', 'scale', 1e0, 'unit', 'a.u.');
format.param.log_var = struct('spec', '%.3f', 'scale', 1e0, 'unit', 'a.u.');
format.param.lin_var = struct('spec', '%.3f', 'scale', 1e0, 'unit', 'a.u.');

end

function [err_mat, wgt_mat] = get_fct_err(param, n_pts)

lin_opt = param.lin_opt;
lin_var = param.lin_var;
log_opt = param.log_opt;
log_var = param.log_var;

% get the error for the log variable
err_lin_mat = (log_var-log_opt).^2;

% get the error for the lin variable
err_log_mat = sum((lin_var-lin_opt).^2, 1);

% add the errors
err_mat = err_lin_mat+err_log_mat;

% single 
wgt_mat = ones(1, n_pts);

end
