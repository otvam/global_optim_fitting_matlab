function [param, optim] = get_solver(var_opt, var_fix, var_err, fct_err, format, cache, optimizer)
% Fit parameters with respect to a dataset with advanced solvers.
%
%    This function provides a common interface for different solvers
%        - fminunc / fminsearch / fmincon
%        - surrogateopt / particleswarm / ga
%        - combination of several solvers
%
%    Customized error function
%        - custom weights for the dataset points
%        - choice of the error norm
%        - recover from undefined values
%        - vectorized evaluation of the error function
%        - caching of the error function
%
%    Advanced parameters handling
%        - scalar or vector parameters
%        - initial values
%        - bounds (even for unconstrained solvers)
%        - variable transformation
%        - variable normlization
%
%    Advanced monitoring capabilities
%        - error metrics
%        - solver metrics
%        - computational cost / timing
%        - plots
%
%    Warning:
%        All the provided features have a computational cost.
%        Therefore, this library is mostly adapted to error functions
%        with large computational cost and/or vectorized error functions.
%
%    Parameters:
%        var_opt (cell): description of the parameters to be fitted
%        var_fix (cell): description of the parameters with fixed values
%        var_err (struct): description of the used error metric
%        fct_err (handle): error function for determining the parameters
%        format (struct): structure with formatting instructions (name and unit)
%        cache (struct): structure with the cache options
%        optimizer (cell): structure describing the solver parameters
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
fprintf('get var\n')
obj_var = SolverVar(var_opt, var_fix, var_err);

% object interfacing the different solvers
fprintf('get interface\n')
obj_run = SolverRun(obj_var, fct_err, format, cache);

% calling the different solvers (using the results as initial values)
disp('run solvers')
[n_pts, param] = obj_var.get_init();
for i=1:length(optimizer)
    [n_pts, param, optim{i}] = obj_run.get_run(n_pts, param, optimizer{i});
end

% solution should be a single parameter combination
assert(n_pts==1, 'invalid solution')

end
