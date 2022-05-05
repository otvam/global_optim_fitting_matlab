classdef SolverRun < handle
    % Class calling the different solvers.
    %
    %    Abstraction layers for the different solvers.
    %    Manage the cache object.
    %    Manage the variable object.
    %    Create and manage the log object. 
    %
    %    Thomas Guillod.
    %    2021-2022 - BSD License.
    
    %% properties
    properties (SetAccess = private, GetAccess = private)
        obj_var % object managing the variables
        obj_cache % object managing the cache for the error function
        format % structure with formatting instructions (name and unit)
    end
    
    %% public
    methods (Access = public)
        function self = SolverRun(obj_var, obj_cache, format)
            % Constructor.
            
            % set data
            self.obj_var = obj_var;
            self.obj_cache = obj_cache;
            self.format = format;
            
            % disable warning if parallel is not required
            warning('off', 'optimlib:commonMsgs:NoPCTLicense');
        end
        
        function [x_scale, optim] = get_run(self, x0_scale, optimizer)
            % Call a specific solver with given initial values.
                        
            % extract
            solver_type = optimizer.solver_type;
            log_iter = optimizer.log_iter;
            log_final = optimizer.log_final;
            clamp_bnd = optimizer.clamp_bnd;
            recover_val = optimizer.recover_val;
            options = optimizer.options;
            
            % get the error function from the cache
            self.obj_cache.get_clear();
            fct_err = @(x_scale) self.obj_cache.get_eval(x_scale);

            % functions for scaling the variables
            fct_unclamp = @(x0_scale) self.obj_var.get_unclamp(x0_scale, clamp_bnd);
            fct_clamp = @(x_clamp) self.obj_var.get_clamp(x_clamp, clamp_bnd);
            fct_param = @(x_scale) self.obj_var.get_param(x_scale);
            
            % get the error function
            fct_sol = @(x_unclamp) SolverRun.get_sol(x_unclamp, fct_err, fct_clamp, recover_val);
            
            % data structure for the logging
            data_optim.fct_err = fct_err;
            data_optim.fct_clamp = fct_clamp;
            data_optim.fct_param = fct_param;
            data_optim.n_var = size(x0_scale, 2);
            
            % create the logging object
            obj_log = SolverLog(solver_type, log_iter, log_final, data_optim, self.format);

            % call the solver
            [x_scale, optim] = SolverRun.get_solver(fct_sol, fct_unclamp, fct_clamp, x0_scale, options, solver_type, obj_log);
        end
    end
    
    %% private static api
    methods (Static, Access = private)
        function [x_scale, optim] = get_solver(fct_sol, fct_unclamp, fct_clamp, x0_scale, options, solver_type, obj_log)
            % Call a solver and manage the logging.
            
            % logging function to be called after each solver iteration
            fct_iter = @(x_unclamp, err, n_iter, n_eval, msg) obj_log.get_iter(x_unclamp, err, n_iter, n_eval, msg);
            
            % logging function to be called after the final solver iteration
            fct_final = @(x_unclamp, err, n_iter, n_eval, msg, is_valid) obj_log.get_final(x_unclamp, err, n_iter, n_eval, msg, is_valid);
            
            % transform bounded variables into unconstrained variables with sine transformation
            [x0_unclamp, lb_unclamp, ub_unclamp] = fct_unclamp(x0_scale);
            
            % run the solver
            x_unclamp = SolverList.get_solver(fct_sol, fct_iter, fct_final, x0_unclamp, lb_unclamp, ub_unclamp, options, solver_type);
            
            % transform unconstrained variables into bounded variables with sine transformation
            x_scale = fct_clamp(x_unclamp);
            
            % get the logging data
            optim = obj_log.get_optim();
        end
        
        function err = get_sol(x_unclamp, fct_err, fct_clamp, recover_val)
            % Error function that will be called by the different solvers.
            
            % transform unconstrained variables into bounded variables with sine transformation
            x_scale = fct_clamp(x_unclamp);
            
            % call the error function
            err = fct_err(x_scale);
                                                                                                
            % replace bad values by the specified value
            idx = isfinite(err)==false;
            err(idx) = recover_val;
        end
    end
end