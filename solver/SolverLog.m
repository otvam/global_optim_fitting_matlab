classdef SolverLog < handle
    % Class for displaying and logging the solver progress / results.
    %
    %    Log and display the solver progress after each iteration.
    %    Log, display, and plot the solver results after the final iteration.
    %
    %    Thomas Guillod.
    %    2021-2022 - BSD License.
    
    %% properties
    properties (SetAccess = private, GetAccess = private)
        solver_type % name of the used solver
        iter % log (or not) the solver iterations
        final % log (or not) the solver final results
        fct_unscale % function the extract the parameter structure from a raw matrix
        fct_sol % function returning the error metrics
        
        t_start % timestamp set at the initialization
        t_last % timestamp of the last iteration
        optim_iter % cell containing the logged solver progress data
        optim_final % cell containing the final iteration data
        obj_disp % object managing the plots and console display
    end
    
    %% public
    methods (Access = public)
        function self = SolverLog(solver_type, log, fct_unscale, fct_sol, format)
            % Constructor.
            
            % set data
            self.solver_type = solver_type;
            self.iter = log.iter;
            self.final = log.final;
            self.fct_unscale = fct_unscale;
            self.fct_sol = fct_sol;
            
            % init the timing
            self.t_start = datetime('now');
            self.t_last = datetime('now');
            
            % init the logging data
            self.optim_iter = {};
            self.optim_final = struct();
            
            % object managing the plots and console display
            self.obj_disp = SolverDisp(solver_type, format);
        end
        
        function get_iter(self, x_scale, err, n_iter, n_eval, msg, is_valid)
            % Log and display the solver progress after an iteration.
            
            if self.iter.log==true
                % get the total solver timing
                [t_solver, t_iter] = get_time(self, n_iter, false);
                
                % get message
                msg = sprintf('intermediate results\nstate: %s', msg);
                
                % get and add the logging data
                optim_tmp = self.get_log(x_scale, err, is_valid, n_iter, n_eval, msg, t_solver, t_iter);
                self.optim_iter{end+1} = optim_tmp;
                
                % display and plot the data
                if self.iter.display==true
                    self.obj_disp.get_disp('iter', n_iter, self.optim_iter{end})
                end
                if self.iter.plot==true
                    self.obj_disp.get_plot_single('iter', n_iter, optim_tmp)
                    self.obj_disp.get_plot_all('iter', n_iter, self.optim_iter)
                    drawnow();
                end
            end
        end
        
        function get_final(self, x_scale, err, n_iter, n_eval, msg, is_valid)
            % Log, display, and plot the solver results after the final iteration.
            
            if self.final.log==true
                % get the total solver timing
                [t_solver, t_iter] = get_time(self, n_iter, true);
                
                % get message
                msg = sprintf('final results\n%s', msg);
                
                % get and add the logging data
                optim_tmp = self.get_log(x_scale, err, is_valid, n_iter, n_eval, msg, t_solver, t_iter);
                self.optim_final = optim_tmp;
                
                % display and plot the data
                if self.final.display==true
                    self.obj_disp.get_disp('final', n_iter, self.optim_iter{end})
                end
                if self.final.plot==true
                    self.obj_disp.get_plot_single('final', n_iter, optim_tmp)
                    self.obj_disp.get_plot_all('final', n_iter, self.optim_iter)
                    drawnow();
                end
            end
        end
        
        function optim = get_optim(self)
            % Get the logged data.
            
            optim.solver_type = self.solver_type;
            optim.iter = self.optim_iter;
            optim.final = self.optim_final;
        end
    end
    
    %% private api
    methods( Access = private)
        function [t_solver, t_iter] = get_time(self, n_iter, is_final)
            % Compute and update the solver timing (total time and iteration time).
            
            % get current time
            t_now = datetime('now');
            
            % get total solver time
            t_solver = t_now-self.t_start;
            
            % the time of the last iteration or average time per iteration
            if is_final==true
                t_iter = t_solver./(n_iter+1);
            else
                t_iter = t_now-self.t_last;
            end
            
            % update the iteration timestamp
            self.t_last = t_now;
        end
        
        function optim = get_log(self, x_scale, err, is_valid, n_iter, n_eval, msg, t_solver, t_iter)
            % Parse and assign the logging data for a specific iteration.
            
            % check if there is a valid solution
            has_solution = (isempty(x_scale)==false)&&(isempty(err)==false)&&any(isfinite(err));
            
            % parse the best solution (if any)
            if has_solution==true
                % best the best parameter combination
                [~, idx_best] = min(err);
                x_scale = x_scale(:,idx_best);
                
                % get the error metrics
                [err_best, err_vec, wgt_vec] = self.fct_sol(x_scale);
                
                % invalid solution can be there due to cache tolerance
                has_solution = isfinite(err_best)&&all(isfinite(err_best))&&all(isfinite(wgt_vec));
                
                % parse data if a solution exists
                if has_solution==true
                    % extract the parameter structure from a raw data
                    [n_pts, param, bnd, is_bound] = self.fct_unscale(x_scale);
                    assert(n_pts==1, 'invalid size: solution')
                    
                    % process the error metrics
                    err_fom = SolverLog.get_err_fom(err_best, err_vec, wgt_vec);
                end
            end
            
            % create a empty solution if required
            if has_solution==false
                param = [];
                bnd = [];
                err_fom = [];
                err_best = NaN;
                is_valid = false;
                is_bound = false;
            end
            
            % get the solver figures of merit
            sol_fom = SolverLog.get_sol_fom(err, err_best, is_valid, is_bound, n_iter, n_eval, msg, t_solver, t_iter);
            
            % assign base data
            optim.has_solution = has_solution;
            optim.sol_fom = sol_fom;
            optim.param = param;
            optim.bnd = bnd;
            optim.err_fom = err_fom;
        end
    end
    
    %% private static api
    methods(Static, Access = private)
        function sol_fom = get_sol_fom(err, err_best, is_valid, is_bound, n_iter, n_eval, msg, t_solver, t_iter)
            % Process the solver data and assign the results to a struct.
            
            % check population
            pop_valid = isfinite(err);
            n_pop_all = length(pop_valid);
            n_pop_valid = nnz(pop_valid);
            
            % parse the solver message
            msg = strtrim(msg);
            msg = splitlines(msg);
            msg = strtrim(msg);
            
            % get the solver figures of merit
            sol_fom.msg = msg;
            sol_fom.err_best = err_best;
            sol_fom.is_bound = is_bound;
            sol_fom.is_valid = is_valid;
            sol_fom.n_iter = n_iter;
            sol_fom.n_eval = n_eval;
            sol_fom.n_pop_all = n_pop_all;
            sol_fom.n_pop_valid = n_pop_valid;
            sol_fom.t_solver = t_solver;
            sol_fom.t_iter = t_iter;
        end
        
        function err_fom = get_err_fom(err_best, err_vec, wgt_vec)
            % Process the error metrics and assign the results to a struct.
            
            % get the error for different types of norms
            err_fom.n_set = (length(err_vec)+length(wgt_vec))./2;
            err_fom.wgt_sum = sum(wgt_vec);
            err_fom.wgt_avg = sum(wgt_vec)./length(wgt_vec);
            
            if err_fom.n_set>1
                % get error with simple metrics
                err_fom.avg = SolverUtils.get_error(err_vec, wgt_vec, 'avg');
                err_fom.min = SolverUtils.get_error(err_vec, wgt_vec, 'min');
                err_fom.max = SolverUtils.get_error(err_vec, wgt_vec, 'max');
                
                % get error for different norm
                norm = [1 2 4 6 8 10 12 Inf];
                for i=1:length(norm)
                    err_fom.norm_err(i) = SolverUtils.get_norm(err_vec, wgt_vec, norm(i));
                    err_fom.norm_val(i) = norm(i);
                end
                
                % get error for different percentiles
                percentile = [0.01 0.05 0.1 0.5 0.9 0.95 0.99];
                for i=1:length(percentile)
                    err_fom.percentile_err(i) = SolverUtils.get_percentile(err_vec, wgt_vec, percentile(i));
                    err_fom.percentile_val(i) = percentile(i);
                end
            end
            
            % assign raw data
            err_fom.err_best = err_best;
            err_fom.err_vec = err_vec;
            err_fom.wgt_vec = wgt_vec;
        end
    end
end