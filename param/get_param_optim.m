function [var_opt, var_fix, var_err, fct_err, format] = get_param_optim()
% Definition of a dummy optimization problem (variables and objective).
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

% definition of the variables optimized by the solver
%    - name (str): name of the variable
%    - x0 (matrix): initial values of the variable (row: size of the variable / col: number of initial values)
%    - idx (vector): index of the variable to be assigned (empty for all)
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
var_opt{end+1} = struct('name', 'log_var', 'x0', 1.0, 'idx', [], 'lb', 0.01, 'ub', 100.0, 'trf', 'log', 'norm', true, 'tol_bnd', 0.05);
var_opt{end+1} = struct('name', 'lin_var', 'x0', [2.0 ; 2.0], 'idx', [1 ; 3], 'lb', 1.0, 'ub', 3.0, 'trf', 'lin', 'norm', true, 'tol_bnd', 0.05);
var_opt{end+1} = struct('name', 'lin_var', 'x0', [2.0 ; 2.0], 'idx', [2 ; 4], 'lb', 0.0, 'ub', 2.0, 'trf', 'lin', 'norm', true, 'tol_bnd', 0.05);

% definition of the variables with constant values (not optimized)
%    - name (str): name of the variable
%    - x0 (column vector): value of the variable (can be a column vector)
%    - idx (vector): index of the variable to be assigned (empty for all)
var_fix = {};
var_fix{end+1} = struct('name', 'log_opt', 'x0', 15.0, 'idx', []);
var_fix{end+1} = struct('name', 'lin_opt', 'x0', [1.5 ; 2.5], 'idx', [1 ; 3]);
var_fix{end+1} = struct('name', 'lin_opt', 'x0', [0.5 ; 1.0], 'idx', [2 ; 4]);

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
var_err = struct('type', 'error', 'arg', 'avg');

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
format.err = struct('spec', '%.3g', 'scale', 1e0, 'unit', 'a.u.', 'xscale', 'lin', 'yscale', 'log');

% formatting instruction for displaying the variables (variables can be omitted)
%    - spec (str): fprintf format specification
%    - scale (str): scaling factor
%    - unit (str): unit of the variable
format.param.log_opt = struct('spec', '%.3f', 'scale', 1e0, 'unit', 'a.u.');
format.param.lin_opt = struct('spec', '%.3f', 'scale', 1e0, 'unit', 'a.u.');
format.param.log_var = struct('spec', '%.3f', 'scale', 1e0, 'unit', 'a.u.');
format.param.lin_var = struct('spec', '%.3f', 'scale', 1e0, 'unit', 'a.u.');

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
%    This example is using the error function as a objective function:
%        - the error function return a scalar value (n_set = 1)
%        - the error function will be minized by the solver
%
%    Parameters:
%        param (struct): parameters combination to be evaluated
%        n_pts (int): number of parameter combinations (vectorized call)
%
%    Returns:
%        n_set (int): number of points in the error / weight
%        err_mat (matrix): error matrix
%        wgt_mat (matrix): weight matrix

% extract the constant
log_opt = param.log_opt;
lin_opt = param.lin_opt;

% extract the variables
log_var = param.log_var;
lin_var = param.lin_var;

% get the error for the log variable (dummy function)
err_log_mat = (log10(log_var)-log10(log_opt)).^2;

% get the error for the lin variable (dummy function)
err_lin_mat = sum((lin_var-lin_opt).^2, 1);

% add the errors (dummy function)
err_mat = err_lin_mat+err_log_mat;

% this is a scalar optimization problem
n_set = 1;

% set uniform weights
wgt_mat = ones(n_set, n_pts);

end
