function [var_opt, var_fix, var_err, fct_err, format] = get_param_fitting()

%% get the dataset
[n_set, val, wgt, x, y] = get_dataset();

%% variables
var_fix = {};
var_fix{end+1} = struct('name', 'n_set', 'x0', n_set);
var_fix{end+1} = struct('name', 'wgt', 'x0', wgt);
var_fix{end+1} = struct('name', 'x', 'x0', x);
var_fix{end+1} = struct('name', 'y', 'x0', y);
var_fix{end+1} = struct('name', 'val', 'x0', val);

var_opt = {};
var_opt{end+1} = struct('name', 'k', 'x0', [0.5, 0.6], 'lb', 0.1, 'ub', 100.0, 'trf', 'log', 'norm', true, 'tol_bnd', 0.05);
var_opt{end+1} = struct('name', 'dx', 'x0', [1.5, 2.5], 'lb', 1.0, 'ub', 3.0, 'trf', 'lin', 'norm', true, 'tol_bnd', 0.05);
var_opt{end+1} = struct('name', 'dy', 'x0', [1.5, 2.5], 'lb', 1.0, 'ub', 3.0, 'trf', 'lin', 'norm', true, 'tol_bnd', 0.05);

var_err = struct('type', 'norm', 'arg', 8);

%% error function
fct_err = @(param, n_pts) get_fct_err(param, n_pts);

%% format
format.err = struct('spec', '%.3f', 'scale', 1e2, 'unit', '%');
format.param.n_set = struct('spec', '%d', 'scale', 1e0, 'unit', '#');
format.param.k = struct('spec', '%.3g', 'scale', 1e0, 'unit', 'a.u.');
format.param.dx = struct('spec', '%.3g', 'scale', 1e0, 'unit', 'a.u.');
format.param.dy = struct('spec', '%.3g', 'scale', 1e0, 'unit', 'a.u.');

end

function [n_set, err_mat, wgt_mat] = get_fct_err(param, n_pts)

% extract the variables
n_set = param.n_set;
x = param.x;
y = param.y;
val = param.val;
wgt = param.wgt;
k = param.k;
dx = param.dx;
dy = param.dy;

% reshape with the number of parameter combinations
val_mat = repmat(val, 1, n_pts);
wgt_mat = repmat(wgt, 1, n_pts);

% evaluate the model
val_model_mat = k.*(x.^dx).*(y.^dy);

% get the error between the model and the dataset
err_mat = (val_model_mat-val_mat)./val_mat;

end

function [n_set, val, wgt, x, y] = get_dataset()

% get the points composing the dataset
x_vec = logspace(log10(1e3), log10(100e3), 25);
y_vec = logspace(log10(1e-2), log10(1e0), 25);
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

% dataset size
n_set = length(wgt);

% the parameters used to generate the value are noisy
k = 15.0.*(1+0.1.*randn(n_set, 1));
dx = 1.5;
dy = 2.5;

% generate the value of the dataset
val = k.*(x.^dx).*(y.^dy);

end
