function [var_opt, var_fix, var_err, fct_err, format] = get_param_optim()

var_fix = {};
var_fix{end+1} = struct('name', 'cst_scalar', 'x0', 0.5, 'idx', 1);
var_fix{end+1} = struct('name', 'cst_vector', 'x0', 1.5, 'idx', 1);
var_fix{end+1} = struct('name', 'cst_vector', 'x0', 2.5, 'idx', 2);

var_opt = {};
var_opt{end+1} = struct('name', 'var_scalar', 'x0', 0.5, 'lb', 0.1, 'ub', 1.0, 'tol_bnd', 0.05, 'trf', 'log', 'norm', true, 'idx', 1);
var_opt{end+1} = struct('name', 'var_vector', 'x0', 2.0, 'lb', 1.0, 'ub', 3.0, 'tol_bnd', 0.05, 'trf', 'lin', 'norm', true, 'idx', 1);
var_opt{end+1} = struct('name', 'var_vector', 'x0', 2.0, 'lb', 1.0, 'ub', 3.0, 'tol_bnd', 0.05, 'trf', 'lin', 'norm', true, 'idx', 2);

var_err = struct('type', 'norm', 'arg', 8);

fct_err = @(param, n) get_fct_err(param, n);

format.err = struct('spec', '%.3f', 'scale', 1e2, 'unit', '%');
format.param.cst_scalar = struct('spec', '%.3g', 'scale', 1e0, 'unit', 'a.u.');
format.param.cst_vector = struct('spec', '%.3g', 'scale', 1e0, 'unit', 'a.u.');
format.param.var_scalar = struct('spec', '%.3g', 'scale', 1e0, 'unit', 'a.u.');
format.param.var_vector = struct('spec', '%.3g', 'scale', 1e0, 'unit', 'a.u.');

end

function [err_mat, wgt_mat] = get_fct_err(param, n)

cst_scalar = param.cst_scalar;
cst_vector = param.cst_vector;
var_scalar = param.var_scalar;
var_vector = param.var_vector;

% compute the value with a dummy model
err_mat = 0;
err_mat = err_mat+(var_scalar-cst_scalar).^2;
err_mat = err_mat+(cst_vector(1,:)-var_vector(1,:)).^2;
err_mat = err_mat+(cst_vector(2,:)-var_vector(2,:)).^2;

wgt_mat = ones(1, n);

end
