function [var_opt, var_fix, fct_err, format] = get_param_problem()

var_fix = {};
var_fix{end+1} = struct('name', 'cst_scalar', 'x', 0.5, 'idx', 1);
var_fix{end+1} = struct('name', 'cst_vector', 'x', 1.0, 'idx', 1);
var_fix{end+1} = struct('name', 'cst_vector', 'x', 1.5, 'idx', 2);

var_opt = {};
var_opt{end+1} = struct('name', 'var_scalar', 'x0', 0.5, 'lb', 0.1, 'ub', 1.0, 'scale', 'log', 'norm', true, 'idx', 1);
var_opt{end+1} = struct('name', 'var_vector', 'x0', 2.0, 'lb', 1.0, 'ub', 3.0, 'scale', 'lin', 'norm', true, 'idx', 1);
var_opt{end+1} = struct('name', 'var_vector', 'x0', 2.0, 'lb', 1.0, 'ub', 3.0, 'scale', 'lin', 'norm', true, 'idx', 2);

cst_scalar = 0.5;
cst_vector = [1.0 ; 1.5];
var_scalar = 0.5;
var_vector = [1.5 ; 2.5];
[val, wgt] = get_model(cst_scalar, cst_vector, var_scalar, var_vector, true);

fct_err = @(param, n) get_fct_err(param, n, val, wgt);

format.err = struct('spec', '%.3f', 'scale', 1e2, 'unit', '%');
format.param.cst_scalar = struct('spec', '%.3g', 'scale', 1e0, 'unit', 'a.u.');
format.param.cst_vector = struct('spec', '%.3g', 'scale', 1e0, 'unit', 'a.u.');
format.param.var_scalar = struct('spec', '%.3g', 'scale', 1e0, 'unit', 'a.u.');
format.param.var_vector = struct('spec', '%.3g', 'scale', 1e0, 'unit', 'a.u.');

end

function [err_mat, wgt_mat] = get_fct_err(param, n, val, wgt)

cst_scalar = param.cst_scalar;
cst_vector = param.cst_vector;
var_scalar = param.var_scalar;
var_vector = param.var_vector;

val_model_mat = get_model(cst_scalar, cst_vector, var_scalar, var_vector, false);

val_mat = repmat(val, 1, n);
wgt_mat = repmat(wgt, 1, n);
err_mat = (val_model_mat-val_mat)./val_mat;

end

function [val, wgt] = get_model(cst_scalar, cst_vector, var_scalar, var_vector, add_noise)

% get the points composing the dataset
x_vec = linspace(1, 5, 10);
y_vec = linspace(1, 5, 10);
[x_mat, y_mat] = ndgrid(x_vec, y_vec);

% get the weigts, double the weight for the edges
wgt_mat = ones(length(x_vec), length(y_vec));
wgt_mat(1,:) = 2;
wgt_mat(end,:) = 2;
wgt_mat(:,1) = 2;
wgt_mat(:,end) = 2;

% flatten the dataset
x = x_mat(:);
y = y_mat(:);
wgt = wgt_mat(:);

% add noise to the dataset points
if add_noise==true
    noise = 0.005;
    x = x.*(1+noise.*(2.*rand(size(x))-1));
    y = y.*(1+noise.*(2.*rand(size(x))-1));
end

% compute the value with a dummy model
val = 0;
val = val+var_scalar.*(cst_scalar+x.^2+y.^2);
val = val+abs(x+cst_vector(1,:)).^var_vector(1,:);
val = val+abs(y+cst_vector(2,:)).^var_vector(2,:);

end
