function [var_opt, var_fix, fct_vec] = get_param_problem()

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
[val, wgt] = get_model_vec(cst_scalar, cst_vector, var_scalar, var_vector, true);

fct_vec = @(param, n) get_fct_vec(param, n, val, wgt);

end

function [err_vec, wgt_vec] = get_fct_vec(param, n, val, wgt)

cst_scalar = param.cst_scalar;
cst_vector = param.cst_vector;
var_scalar = param.var_scalar;
var_vector = param.var_vector;

val_model_vec = get_model_vec(cst_scalar, cst_vector, var_scalar, var_vector, false);

val_vec = repmat(val, 1, n);
wgt_vec = repmat(wgt, 1, n);
err_vec = (val_model_vec-val_vec)./val_vec;

end

function [val, wgt] = get_model_vec(cst_scalar, cst_vector, var_scalar, var_vector, add_noise)

% get grid
x_vec = linspace(1, 5, 10);
y_vec = linspace(1, 5, 10);
noise = 0.1;

% get point and weights
[x_mat, y_mat] = ndgrid(x_vec, y_vec);
wgt_mat = ones(length(x_vec), length(y_vec));
wgt_mat(1,:) = 2;
wgt_mat(end,:) = 2;
wgt_mat(:,1) = 2;
wgt_mat(:,end) = 2;

% flatten
x = x_mat(:);
y = y_mat(:);
wgt = wgt_mat(:);

% add noise
if add_noise==true
    x = x.*(1+noise.*(2.*rand(size(x))-1));
    y = y.*(1+noise.*(2.*rand(size(x))-1));
end

val = 0;
val = val+var_scalar.*(cst_scalar+x.^2+y.^2);
val = val+abs(x+cst_vector(1,:)).^var_vector(1,:);
val = val+abs(y+cst_vector(2,:)).^var_vector(2,:);

end
