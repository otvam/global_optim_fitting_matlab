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
        fct_err % error function for determining the parameters
        format % structure with formatting instructions (name and unit)
    end
    
    %% public
    methods (Access = public)
        function self = SolverRun(obj_var, fct_err, format)
            % Constructor.
            
            % set data
            self.obj_var = obj_var;
            self.fct_err = fct_err;
            self.format = format;
            
            % disable warning if parallel is not required
            warning('off', 'optimlib:commonMsgs:NoPCTLicense');
        end
        
        function [n_pts, param, optim] = get_run(self, n_pts, param, optimizer)
            % Call a specific solver with given initial values.
            
            % extract
            solver_type = optimizer.solver_type;
            log_iter = optimizer.log_iter;
            log_final = optimizer.log_final;
            clamp_bnd = optimizer.clamp_bnd;
            recover_val = optimizer.recover_val;
            options = optimizer.options;
            
            % get the error function from the cache
            %             self.obj_cache.get_clear();
            %             fct_err = @(x_scale) self.obj_cache.get_eval(x_scale);
            
            % get variables scaling
            [n_var, x_scale, lb_scale, ub_scale] = self.obj_var.get_scale(n_pts, param, clamp_bnd);
            fct_unscale = @(x_scale) self.obj_var.get_unscale(x_scale, clamp_bnd);
            fct_scale_err = @(err_mat, wgt_mat) self.obj_var.get_scale_err(err_mat, wgt_mat);
            
            % get the error function
            fct_sol = @(x_scale) SolverRun.get_sol(x_scale, self.fct_err, fct_unscale, fct_scale_err);
            fct_recover = @(x_scale) SolverRun.get_sol_recover(x_scale, fct_sol, recover_val);
                        
            % data structure for the logging
            data_optim.fct_sol = fct_sol;
            data_optim.fct_unscale = fct_unscale;
            data_optim.n_var = n_var;
            
            % create the logging object
            obj_log = SolverLog(solver_type, log_iter, log_final, data_optim, self.format);
            
            % call the solver
            [x_scale, optim] = SolverRun.get_solver(fct_recover, x_scale, lb_scale, ub_scale, options, solver_type, obj_log);
            
            % unscale
            [n_pts, param] = fct_unscale(x_scale);
        end
    end
    
    %% private static api
    methods (Static, Access = private)
        function [x_scale, optim] = get_solver(fct_recover, x_scale, lb_scale, ub_scale, options, solver_type, obj_log)
            % Call a solver and manage the logging.
            
            % logging function to be called after each solver iteration
            fct_iter = @(x_scale, err, n_iter, n_eval, msg) obj_log.get_iter(x_scale.', err.', n_iter, n_eval, msg);
            
            % logging function to be called after the final solver iteration
            fct_final = @(x_scale, err, n_iter, n_eval, msg, is_valid) obj_log.get_final(x_scale.', err.', n_iter, n_eval, msg, is_valid);
            
            % run the solver
            x_scale = x_scale.';
            lb_scale = lb_scale.';
            ub_scale = ub_scale.';
            x_scale = SolverList.get_solver(fct_recover, fct_iter, fct_final, x_scale, lb_scale, ub_scale, options, solver_type);
            x_scale = x_scale.';
            
            % get the logging data
            optim = obj_log.get_optim();
        end
                
        function [err_mat, wgt_mat] = get_cache(x_scale, fct_err, fct_unscale)
                        
            % unscale variables
            [n_pts, param] = fct_unscale(x_scale);

            % call the error function
            [err_mat, wgt_mat] = fct_err(param, n_pts);
            
        end
        
        function [err_best, n_set, err_mat, wgt_mat] = get_sol(x_scale, fct_err, fct_unscale, fct_scale_err)
            % Error function that will be called by the different solvers.
            
            % evaluate the function
            if isempty(x_scale)
                err_best = [];
                n_set = [];
                err_mat = [];
                wgt_mat = [];
            else
                [err_mat, wgt_mat] = SolverRun.get_cache(x_scale, fct_err, fct_unscale);
                [err_best, n_set] = fct_scale_err(err_mat, wgt_mat);
            end
        end
        
        function err_best = get_sol_recover(x_scale, fct_sol, recover_val)
            
            % reshape
            x_scale = x_scale.';
            
            % get error
            err_best = fct_sol(x_scale);
            
            % replace invalid values
            idx = isfinite(err_best)==false;
            err_best(idx) = recover_val;
            
            % reshape
            err_best = err_best.';
        end
    end
end