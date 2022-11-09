function [var_opt, var_fix, var_err, fct_err, format] = get_param_fitting()
% Definition of a dummy fitting problem (variables and objective).
%
%    Returns:
%        var_opt (cell): description of the parameters to be fitted
%        var_fix (cell): description of the parameters with fixed values
%        var_err (struct): description of the used error metric
%        fct_err (handle): error function for determining the parameters
%        format (struct): structure with formatting instructions (name and unit)
%
%    Thomas Guillod.
%    2021-2022 - BSD License.

% get the dataset that will be used for the fitting
[n_set, val, wgt, x, y] = get_dataset();

% definition of the variables optimized by the solver
%    - name (str): name of the variable
%    - idx (vector): index of the variable to be assigned (empty for all)
%    - x0 (matrix): initial values of the variable (row: size of the variable / col: number of initial values)
%    - lb (scalar): lower bound of the variable (infinite value are tolerated)
%    - ub (scalar): upper bound of the variable (infinite value are tolerated)
%    - trf (str): apply a variable transformation for the solver to improve conditionning
%        - 'lin': linear variable transformation (do nothing)
%        - 'sqrt': consider the square root of the variable
%        - 'quad': consider the square of the variable
%        - 'log': consider the log of the variable
%    - norm (bool): normalize the variable for the solver to improve conditionning
%        - the variable is scaled between zero and one using the bounds
%        - this normalization is only used if both bounds are finite
%    - tol_bnd (scalar): tolerance for detecting if a variable is close to the bounds
%        - new bounds are selected lb+tol and ub-tol with tol = tol_bnd*(ub-lb)
%        - the solver checks if the variable is between these new bounds
%        - this tolerance is not used if (ub-lb) is not finite
%        - this tolerance is only used for analysing the solution, not for the solver itself
%        - this tolerance is usefull to detect if the solver is stuck near a bound
var_opt = {};
var_opt{end+1} = struct('name', 'k', 'x0', [0.5, 0.6], 'idx', [], 'lb', 0.1, 'ub', 100.0, 'trf', 'log', 'norm', true, 'tol_bnd', 0.05);
var_opt{end+1} = struct('name', 'dx', 'x0', [1.5, 2.5], 'idx', [], 'lb', 1.0, 'ub', 3.0, 'trf', 'lin', 'norm', true, 'tol_bnd', 0.05);
var_opt{end+1} = struct('name', 'dy', 'x0', [1.5, 2.5], 'idx', [], 'lb', 1.0, 'ub', 3.0, 'trf', 'lin', 'norm', true, 'tol_bnd', 0.05);

% definition of the variables with constant values (not optimized)
%    - name (str): name of the variable
%    - x0 (matrix): value of the variable (can be a scalar, a vector, or a matrix)
var_fix = {};
var_fix{end+1} = struct('name', 'n_set', 'x0', n_set);
var_fix{end+1} = struct('name', 'wgt', 'x0', wgt);
var_fix{end+1} = struct('name', 'x', 'x0', x);
var_fix{end+1} = struct('name', 'y', 'x0', y);
var_fix{end+1} = struct('name', 'val', 'x0', val);

% definition of the considered error metric
%    - the error function is returning an error vector and a weight vector
%    - the error vector and the weight vector are combined into a scalar error metric
%    - this scalar error metric is used as the objective function for the solver
%    - the following combinations of 'type' and 'arg' are available
%        - 'error' / 'avg': weighted average
%        - 'error' / 'min': minimum value
%        - 'error' / 'max': maximum value
%        - 'norm' / p:  weighted p-norm (p is between 1 and Inf)
%        - 'percentile' / p:  weighted p-percentile (p is between 0 and 1)
var_err = struct('type', 'norm', 'arg', 8);

% definition of the error function used by the solver
%    - the error function is returning an error vector and a weight vector
%    - many parameter combinations can be evaluated together
fct_err = @(param, n_pts) get_fct_err(param, n_pts);

% number of indent characters for the log display
format.indent = 0;

% formatting instruction for displaying the error metric
%    - spec (str): fprintf format specification
%    - scale (str): scaling factor
%    - unit (str): unit of the variable
%    - xscale (str): x-scaling of the convergence plot ('lin' or 'log')
%    - yscale (str): x-scaling of the convergence plot ('lin' or 'log')
format.err = struct('spec', '%.3f', 'scale', 1e2, 'unit', '%', 'xscale', 'lin', 'yscale', 'log');

% formatting instruction for displaying the variables (variables can be omitted)
%    - spec (str): fprintf format specification
%    - scale (str): scaling factor
%    - unit (str): unit of the variable
format.param.n_set = struct('spec', '%d', 'scale', 1e0, 'unit', '#');
format.param.k = struct('spec', '%.3g', 'scale', 1e0, 'unit', 'a.u.');
format.param.dx = struct('spec', '%.3g', 'scale', 1e0, 'unit', 'a.u.');
format.param.dy = struct('spec', '%.3g', 'scale', 1e0, 'unit', 'a.u.');

end

function [n_set, err_mat, wgt_mat] = get_fct_err(param, n_pts)
% Error function used by the solver for finding the optimal parameters.
%
%    The function return a error vector and weight vector for each parameter combinations.
%    The error vector and weight vector are column vector.
%    The function can be used with vectorized calls.
%
%    The dimension of the matrices are the following:
%        - n_pts: the number of parameter combinations
%        - n_set: the size of the error and weights
%        - err_mat: error matrix (size: n_set x n_pts) 
%        - wgt_mat: weight matrix (size: n_set x n_pts)
%
%    This example is using the error function as a fitting function.
%        - the error function return a vector value (n_set > 1)
%        - each point in the vector and weight vectors represent a dataset point
%
%    Parameters:
%        param (struct): parameters combination to be evaluated
%        n_pts (int): number of parameter combinations (vectorized call)
%
%    Returns:
%        n_set (int): number of points in the error / weight
%        err_mat (matrix): error matrix
%        wgt_mat (matrix): weight matrix

% extract the dataset
n_set = param.n_set;
x = param.x;
y = param.y;
val = param.val;
wgt = param.wgt;

% extract the variables
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
% Create a dummy 2D dataset for the fitting problem.
%
%    Returns:
%        n_set (int): size of the dataset
%        val (vector): value of the dataset points
%        wgt (vector): weight of the dataset points
%        x (vector): x-position of the dataset points
%        y (vector): y-position of the dataset points

% get the points composing the dataset
x_vec = logspace(log10(1e3), log10(100e3), 25);
y_vec = logspace(log10(1e-2), log10(1e0), 25);
[x_mat, y_mat] = ndgrid(x_vec, y_vec);

% get the weigts, increased weight for the edges
wgt_mat = NaN(length(x_vec), length(y_vec));
for i=0:8
    wgt_mat((1+i):(end-i), (1+i):(end-i)) = 9-i;
end

% flatten the dataset
x = x_mat(:);
y = y_mat(:);
wgt = wgt_mat(:);

% dataset size
n_set = length(wgt);

% the parameters used to generate the value are noisy
k = 15.0.*(1+0.01.*randn(n_set, 1));
dx = 1.5;
dy = 2.5;

% generate the value of the dataset
val = k.*(x.^dx).*(y.^dy);

end
