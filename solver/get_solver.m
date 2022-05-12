function [param, optim] = get_solver(var_opt, var_fix, var_err, fct_err, format, cache, optimizer)
% MATLAB Toolbox for Global Fitting/Optimization.
%
%    This MATLAB toolbox can be used for the following problems:
%        - finding global minimum of a function
%        - fitting a function to a dataset
%
%    This toolbox is specially adapted to the following problems:
%        - non-smooth error function
%        - non-convex error function
%        - computationally heavy error function
%        - error function with local minima
%        - error function with many input variables
%
%    This toolbox provides a common interface for different solvers:
%        - gradient: fminunc / fmincon
%        - simplex: fminsearch
%        - surrogate: surrogateopt
%        - evolutionary: particleswarm / ga
%        - the aforementioned solvers can be combined
%
%    Customized error function:
%        - custom weights for the dataset points
%        - choice of the error metric (norm, average, percentile, etc.)
%        - recover from undefined values
%        - vectorized evaluation of the error function
%        - parallel evaluation of the error function
%        - caching of the error function
%
%    Advanced variable handling:
%        - initial values
%        - scalar or vector variables
%        - variable transformation (linear, quadratic, logarithmic, etc.)
%        - variable normalization
%        - constraints (lower and upper bounds)
%        - transformation of constrained variables in unconstrained variables
%
%    Advanced monitoring capabilities:
%        - compute various error metrics
%        - compute solver figures of merit
%        - plot/display the solver progress
%        - plot/display the final results
%
%    Warning:
%        All the provided features have a computational cost.
%        Therefore, this library is mostly adapted to time-consuming error functions.
%        For simple error functions, the overhead is non-negligible.
%
%    Parameters:
%        var_opt (cell): description of the parameters to be fitted
%        var_fix (cell): description of the parameters with fixed values
%        var_err (struct): description of the used error metric
%        fct_err (handle): error function for determining the parameters
%        format (struct): structure with formatting instructions (name and unit)
%        cache (struct): structure with the cache options
%        optimizer (cell): cell describing the different solvers to be used
%
%    Returns:
%        param (struct): fitted parameters
%        optim (cell): solver metrics
%
%    Thomas Guillod.
%    2021-2022 - BSD License.

% disable warning if parallel is not required
warning('off', 'optimlib:commonMsgs:NoPCTLicense');

% object managing the variables (name, bounds, transformation, normalization, etc.)
obj_var = SolverVar(var_opt, var_fix, var_err);

% object interfacing the different solvers
obj_run = SolverRun(obj_var, fct_err, format, cache);

% calling the different solvers (using the results as initial values)
[n_pts, param] = obj_var.get_init();
for i=1:length(optimizer)
    [n_pts, param, optim{i}] = obj_run.get_run(n_pts, param, optimizer{i});
end

% solution should be a single parameter combination
assert(n_pts==1, 'invalid solution')

end
