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
        format % structure with formatting instructions (name and unit)
        
        t_start % timestamp set at the initialization
        t_last % timestamp of the last iteration
        optim % cell containing the logged data
        handle_single % figure containing the info on a iteration
        handle_all % figure containing the info on all iterations
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
            self.format = format;
            
            % init the timing and logging data
            self.t_start = datetime('now');
            self.t_last = datetime('now');
            self.optim = {};
            
            % get figure handles
            if (self.iter.plot==true)||(self.final.plot==true)
                self.handle_single = SolverLog.get_figure([self.solver_type ' / single'], 2);
                self.handle_all = SolverLog.get_figure([self.solver_type ' / all'], 1);
            else
                self.handle_single = [];
                self.handle_all = [];
            end
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
                self.optim{end+1} = optim_tmp;
                
                % display the data
                name = sprintf('iter / %s / n = %d', self.solver_type, n_iter);
                
                if self.iter.display==true
                    SolverLog.get_disp(name, self.optim{end}, self.format)
                end
                if self.iter.plot==true
                    SolverLog.get_plot_single(self.handle_single, name, self.optim{end}, self.format)
                    SolverLog.get_plot_all(self.handle_all, name, self.optim, self.format)
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
                self.optim{end+1} = optim_tmp;
                
                % display and plot the data
                name = sprintf('final / %s / n = %d', self.solver_type, n_iter);
                if self.final.display==true
                    SolverLog.get_disp(name, self.optim{end}, self.format)
                end
                if self.final.plot==true
                    SolverLog.get_plot_single(self.handle_single, name, self.optim{end}, self.format)
                    SolverLog.get_plot_all(self.handle_all, name, self.optim, self.format)
                    drawnow();
                end
            end
        end
        
        function optim = get_optim(self)
            % Get the logged data.
            
            optim = self.optim;
        end
    end
    
    %% public static api
    methods(Static, Access = public)
        function get_disp(name, optim, format)
            % Display the logged data for a specific iteration.
            
            % extract
            err = format.err;
            param = format.param;
            indent = format.indent;

            % name of the iteration
            pad = repmat(' ', 1, indent);
            SolverLog.get_print(pad, '%s', name)
            
            % display the solver figures of merit
            SolverLog.get_print(pad, '    fom')
            SolverLog.get_print(pad, '        msg')
            for i=1:length(optim.sol_fom.msg)
                SolverLog.get_print(pad, '            %s', optim.sol_fom.msg{i})
            end
            SolverLog.get_print(pad, '        status')
            SolverLog.get_print(pad, '            is_valid = %s',  mat2str(optim.sol_fom.is_valid))
            SolverLog.get_print(pad, '            is_bound = %s',  mat2str(optim.sol_fom.is_bound))
            SolverLog.get_print(pad, '        count')
            SolverLog.get_print(pad, '            n_iter = %d', optim.sol_fom.n_iter)
            SolverLog.get_print(pad, '            n_eval = %d', optim.sol_fom.n_eval)
            SolverLog.get_print(pad, '        timing')
            SolverLog.get_print(pad, '            t_solver = %s', char(optim.sol_fom.t_solver))
            SolverLog.get_print(pad, '            t_iter = %s', char(optim.sol_fom.t_iter))
            SolverLog.get_print(pad, '        error')
            SolverLog.get_print(pad, '            n_pop_all = %d', optim.sol_fom.n_pop_all)
            SolverLog.get_print(pad, '            n_pop_fail = %d', optim.sol_fom.n_pop_fail)
            SolverLog.get_print(pad, '            err_best = %s', SolverLog.get_format_scalar(optim.sol_fom.err_best, err))
            
            if optim.has_solution==true
                % display the error metrics
                SolverLog.get_print(pad, '    err_fom')
                SolverLog.get_print(pad, '        size')
                SolverLog.get_print(pad, '            n_set = %d', optim.err_fom.n_set)
                SolverLog.get_print(pad, '            err_best = %s', SolverLog.get_format_scalar(optim.sol_fom.err_best, err))
                SolverLog.get_print(pad, '        weight')
                SolverLog.get_print(pad, '            wgt_sum = %.3f', optim.err_fom.wgt_sum)
                SolverLog.get_print(pad, '            wgt_avg = %.3f', optim.err_fom.wgt_avg)
                if optim.err_fom.n_set>1
                    SolverLog.get_print(pad, '        avg')
                    SolverLog.get_print(pad, '            avg = %s', SolverLog.get_format_scalar(optim.err_fom.avg, err))
                    SolverLog.get_print(pad, '            min = %s', SolverLog.get_format_scalar(optim.err_fom.min, err))
                    SolverLog.get_print(pad, '            max = %s', SolverLog.get_format_scalar(optim.err_fom.max, err))
                    SolverLog.get_print(pad, '        norm')
                    for i=1:length(optim.err_fom.norm_val)
                        str_val = sprintf('%d', optim.err_fom.norm_val(i));
                        str_err = SolverLog.get_format_scalar(optim.err_fom.norm_err(i), err);
                        SolverLog.get_print(pad, '            norm / %s = %s', str_val, str_err)
                    end
                    SolverLog.get_print(pad, '        percentile')
                    for i=1:length(optim.err_fom.percentile_val)
                        str_val = sprintf('%.1f %%', 1e2.*optim.err_fom.percentile_val(i));
                        str_err = SolverLog.get_format_scalar(optim.err_fom.percentile_err(i), err);
                        SolverLog.get_print(pad, '            percentile / %s = %s', str_val, str_err)
                    end
                end
                
                % display the current best parameter combinations
                SolverLog.get_print(pad, '    param')
                field = fieldnames(optim.param);
                for i=1:length(field)
                    value = optim.param.(field{i});
                    if isfield(param, field{i})
                        str = SolverLog.get_format_vec(value, param.(field{i}));
                        SolverLog.get_print(pad, '        param.%s = %s', field{i}, str)
                    else
                        SolverLog.get_print(pad, '        param.%s = hidden', field{i})
                    end
                end
                
                % display if the parameters are close to the bounds
                SolverLog.get_print(pad, '    bnd')
                field = fieldnames(optim.bnd);
                for i=1:length(field)
                    value = optim.bnd.(field{i});
                    str = SolverLog.get_format_bnd(value);
                    SolverLog.get_print(pad, '        bnd.%s = %s', field{i}, str)
                end
            end
        end
        
        function get_plot_single(handle, name, optim, format)
            % Plot the logged data for a specific iteration.
            
            if optim.has_solution==true
                % extract format
                scale = format.err.scale;
                unit = format.err.unit;
                
                % extract the error vector and the weighted error vector
                err_vec = optim.err_fom.err_vec;
                wgt_vec = optim.err_fom.wgt_vec;
                n_set = optim.err_fom.n_set;
                
                % plot the error distribution (if it exists)
                if n_set>1
                    if isempty(handle)
                        handle = SolverLog.get_figure(name, 2);
                    end
                    set(handle.fig, 'Visible', 'on')

                    % histogram
                    histogram(handle.ax(1), scale.*err_vec, 'Normalization', 'pdf')
                    grid(handle.ax(1), 'on')
                    xlabel(handle.ax(1), ['errors (' unit ')'], 'interpreter', 'none')
                    ylabel(handle.ax(1), 'probability (#)', 'interpreter', 'none')
                    title(handle.ax(1), sprintf('Histogram / %s', name), 'interpreter', 'none')
                    
                    % scatter plot with error and weights
                    scatter(handle.ax(2), scale.*err_vec, wgt_vec, 50, 'b', 'filled')
                    grid(handle.ax(2), 'on')
                    xlabel(handle.ax(2), ['errors (' unit ')'], 'interpreter', 'none')
                    ylabel(handle.ax(2), 'weights (#)', 'interpreter', 'none')
                    title(handle.ax(2), sprintf('Weights / %s', name), 'interpreter', 'none')
                end
            end
        end
        
        function get_plot_all(handle, name, optim, format)
            % Plot the logged data across all iterations.
            
            % extract format
            scale = format.err.scale;
            unit = format.err.unit;
            xscale = format.err.xscale;
            yscale = format.err.yscale;
            
            % extract the data for all iterations (error metric and timing)
            conv_vec = 1:length(optim);
            for i=1:length(optim)
                err_best = optim{i}.sol_fom.err_best;
                t_solver = optim{i}.sol_fom.t_solver;
                
                err_conv_vec(i) = err_best;
                t_solver_conv_vec(i) = t_solver;
            end
            
            % plot all iterations
            if isempty(handle)
                handle = SolverLog.get_figure(name, 1);
            end
            set(handle.fig, 'Visible', 'on')
            
            % plot the error metric
            plot(handle.ax(1), conv_vec, scale.*err_conv_vec, 'og')
            xlim(handle.ax(1), [min(conv_vec)-1 max(conv_vec)+1])
            grid(handle.ax(1), 'on')
            set(handle.ax(1), 'xscale', xscale)
            set(handle.ax(1), 'yscale', yscale)
            xlabel(handle.ax(1), 'iter (#)', 'interpreter', 'none')
            ylabel(handle.ax(1), ['err (' unit ')'], 'interpreter', 'none')
            title(handle.ax(1), sprintf('Convergence / %s', name), 'interpreter', 'none')
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
                
                % extract the parameter structure from a raw data
                [n_pts, param, bnd, is_bound] = self.fct_unscale(x_scale);
                assert(n_pts==1, 'invalid size: solution')
                
                % get the error metrics
                [err_best, err_vec, wgt_vec] = self.fct_sol(x_scale);
                
                % process the error metrics
                err_fom = SolverLog.get_err_fom(err_best, err_vec, wgt_vec);
            else
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
        function txt = get_format_scalar(val, format)
            % Parse a scalar to a string with scaling and units.
            
            % extract
            spec = format.spec;
            scale = format.scale;
            unit = format.unit;
            
            % parse each elements
            txt = sprintf(spec, scale.*val);
            
            % add unit
            txt = [txt ' ' unit];
        end
        
        function txt = get_format_vec(vec, format)
            % Parse a vector to a string with scaling and units.
            
            % extract
            spec = format.spec;
            scale = format.scale;
            unit = format.unit;
            
            % parse each elements
            for i=1:length(vec)
                txt{i} = sprintf(spec, scale.*vec(i));
            end
            
            % assemble the string
            txt = sprintf('[%s]', strjoin(txt, ' ; '));
            
            % add unit
            txt = [txt ' ' unit];
        end
        
        function txt = get_format_bnd(vec)
            % Parse a boolean vector to string.
            
            % parse each elements
            for i=1:length(vec)
                txt{i} = mat2str(vec(i));
            end
            
            % assemble the string
            txt = sprintf('[%s]', strjoin(txt, ' ; '));
        end
        
        function sol_fom = get_sol_fom(err, err_best, is_valid, is_bound, n_iter, n_eval, msg, t_solver, t_iter)
            % Process the solver data and assign the results to a struct.
            
            % check population
            pop_valid = isfinite(err);
            n_pop_all = length(pop_valid);
            n_pop_fail = nnz(pop_valid==false);
            
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
            sol_fom.n_pop_fail = n_pop_fail;
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
        
        function handle = get_figure(name, n)
            % Create a figure and return axes handles.

            fig = figure('name', name, 'Visible', 'off');
            for i=1:n
                ax(i) = subplot(1, n, i);
            end

            handle.fig = fig;
            handle.ax = ax;
        end
        
        function get_print(pad, name, varargin)
            % Print a line with padding.

            fprintf([pad, name '\n'], varargin{:})
        end
    end
end