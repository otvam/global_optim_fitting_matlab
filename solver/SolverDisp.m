classdef SolverDisp < handle
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
        format % structure with formatting instructions (name and unit)
    end
    
    %% public
    methods (Access = public)
        function self = SolverDisp(solver_type, format)
            % Constructor.
            
            % set data
            self.solver_type = solver_type;
            self.format = format;
        end
                
        function get_disp(self, name, n_iter, optim)
            % Display the logged data for a specific iteration.
            
            % extract
            err = self.format.err;
            param = self.format.param;
            indent = self.format.indent;
            
            % get name
            pad = repmat(' ', 1, indent);
            name = sprintf('%s / %s / n = %d', self.solver_type, name, n_iter);
            
            % name of the iteration
            SolverDisp.get_print(pad, '%s', name)
            
            % display the solver figures of merit
            SolverDisp.get_print(pad, '    fom')
            SolverDisp.get_print(pad, '        msg')
            for i=1:length(optim.sol_fom.msg)
                SolverDisp.get_print(pad, '            %s', optim.sol_fom.msg{i})
            end
            SolverDisp.get_print(pad, '        status')
            SolverDisp.get_print(pad, '            is_valid = %s',  mat2str(optim.sol_fom.is_valid))
            SolverDisp.get_print(pad, '            is_bound = %s',  mat2str(optim.sol_fom.is_bound))
            SolverDisp.get_print(pad, '        count')
            SolverDisp.get_print(pad, '            n_iter = %d', optim.sol_fom.n_iter)
            SolverDisp.get_print(pad, '            n_eval = %d', optim.sol_fom.n_eval)
            SolverDisp.get_print(pad, '        timing')
            SolverDisp.get_print(pad, '            t_solver = %s', char(optim.sol_fom.t_solver))
            SolverDisp.get_print(pad, '            t_iter = %s', char(optim.sol_fom.t_iter))
            SolverDisp.get_print(pad, '            is_timeout = %s',  mat2str(optim.sol_fom.is_timeout))
            SolverDisp.get_print(pad, '        error')
            SolverDisp.get_print(pad, '            n_pop_all = %d', optim.sol_fom.n_pop_all)
            SolverDisp.get_print(pad, '            n_pop_valid = %d', optim.sol_fom.n_pop_valid)
            SolverDisp.get_print(pad, '            err_best = %s', SolverDisp.get_format_scalar(optim.sol_fom.err_best, err))
            
            if optim.has_solution==true
                % display the error metrics
                SolverDisp.get_print(pad, '    err_fom')
                SolverDisp.get_print(pad, '        size')
                SolverDisp.get_print(pad, '            n_set = %d', optim.err_fom.n_set)
                SolverDisp.get_print(pad, '            err_best = %s', SolverDisp.get_format_scalar(optim.sol_fom.err_best, err))
                SolverDisp.get_print(pad, '        weight')
                SolverDisp.get_print(pad, '            wgt_sum = %.3f', optim.err_fom.wgt_sum)
                SolverDisp.get_print(pad, '            wgt_avg = %.3f', optim.err_fom.wgt_avg)
                if optim.err_fom.n_set>1
                    SolverDisp.get_print(pad, '        avg')
                    SolverDisp.get_print(pad, '            avg = %s', SolverDisp.get_format_scalar(optim.err_fom.avg, err))
                    SolverDisp.get_print(pad, '            min = %s', SolverDisp.get_format_scalar(optim.err_fom.min, err))
                    SolverDisp.get_print(pad, '            max = %s', SolverDisp.get_format_scalar(optim.err_fom.max, err))
                    SolverDisp.get_print(pad, '        norm')
                    for i=1:length(optim.err_fom.norm_val)
                        str_val = sprintf('%d', optim.err_fom.norm_val(i));
                        str_err = SolverDisp.get_format_scalar(optim.err_fom.norm_err(i), err);
                        SolverDisp.get_print(pad, '            norm / %s = %s', str_val, str_err)
                    end
                    SolverDisp.get_print(pad, '        percentile')
                    for i=1:length(optim.err_fom.percentile_val)
                        str_val = sprintf('%.1f %%', 1e2.*optim.err_fom.percentile_val(i));
                        str_err = SolverDisp.get_format_scalar(optim.err_fom.percentile_err(i), err);
                        SolverDisp.get_print(pad, '            percentile / %s = %s', str_val, str_err)
                    end
                end
                
                % display the current best parameter combinations
                SolverDisp.get_print(pad, '    param')
                field = fieldnames(optim.param);
                for i=1:length(field)
                    value = optim.param.(field{i});
                    if isfield(param, field{i})
                        str = SolverDisp.get_format_mat(value, param.(field{i}));
                        SolverDisp.get_print(pad, '        param.%s = %s', field{i}, str)
                    else
                        SolverDisp.get_print(pad, '        param.%s = hidden', field{i})
                    end
                end
                
                % display if the parameters are close to the bounds
                SolverDisp.get_print(pad, '    bnd')
                field = fieldnames(optim.bnd);
                for i=1:length(field)
                    value = optim.bnd.(field{i});
                    str = SolverDisp.get_format_bnd(value);
                    SolverDisp.get_print(pad, '        bnd.%s = %s', field{i}, str)
                end
            end
        end
        
        function get_plot_single(self, name, n_iter, optim)
            % Plot the logged data for a specific iteration.
            
            if optim.has_solution==true
                % extract format
                scale = self.format.err.scale;
                unit = self.format.err.unit;
                
                % extract the error vector and the weighted error vector
                err_vec = optim.err_fom.err_vec;
                wgt_vec = optim.err_fom.wgt_vec;
                n_set = optim.err_fom.n_set;
                
                % plot the error distribution (if it exists)
                if n_set>1
                    % get figure and axis handles
                    ax = SolverDisp.get_figure([self.solver_type ' / single'], 2);
                    
                    % get figure name
                    name = sprintf('%s / n = %d', name, n_iter);
                    
                    % histogram
                    histogram(ax(1), scale.*err_vec, 'Normalization', 'pdf')
                    grid(ax(1), 'on')
                    xlabel(ax(1), ['errors (' unit ')'], 'interpreter', 'none')
                    ylabel(ax(1), 'probability (#)', 'interpreter', 'none')
                    title(ax(1), sprintf('Histogram / %s', name), 'interpreter', 'none')
                    
                    % scatter plot with error and weights
                    scatter(ax(2), scale.*err_vec, wgt_vec, 50, scale.*err_vec.*wgt_vec, 'filled')
                    grid(ax(2), 'on')
                    cbar = colorbar(ax(2));
                    cbar.Label.String = ['weights x errors (' unit ')'];
                    cbar.Label.Interpreter = 'none';
                    xlabel(ax(2), ['errors (' unit ')'], 'interpreter', 'none')
                    ylabel(ax(2), 'weights (#)', 'interpreter', 'none')
                    title(ax(2), sprintf('Weights / %s', name), 'interpreter', 'none')
                end
            end
        end
        
        function get_plot_all(self, name, n_iter, optim)
            % Plot the logged data across all iterations.
            
            % extract format
            scale = self.format.err.scale;
            unit = self.format.err.unit;
            xscale = self.format.err.xscale;
            yscale = self.format.err.yscale;
                        
            % extract the data for all iterations
            for i=1:length(optim)
                err_best_vec(i) = optim{i}.sol_fom.err_best;
                n_iter_vec(i) = optim{i}.sol_fom.n_iter;
            end
            
            % get figure and axis handles
            ax = SolverDisp.get_figure([self.solver_type ' / all'], 1);

            % get figure name
            name = sprintf('%s / n = %d', name, n_iter);
            
            % plot the error metric
            plot(ax, n_iter_vec, scale.*err_best_vec, 'or', 'MarkerFaceColor', 'r')
            xlim(ax, [min(n_iter_vec)-1 max(n_iter_vec)+1])
            grid(ax, 'on')
            set(ax, 'xscale', xscale)
            set(ax, 'yscale', yscale)
            xlabel(ax, 'iterations (#)', 'interpreter', 'none')
            ylabel(ax, ['total error (' unit ')'], 'interpreter', 'none')
            title(ax, sprintf('Convergence / %s', name), 'interpreter', 'none')
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
        
        function txt = get_format_mat(mat, format)
            % Parse a matrix to a string with scaling and units.
            
            % extract
            spec = format.spec;
            scale = format.scale;
            unit = format.unit;
            
            % parse each elements
            txt = {};
            for i=1:size(mat, 1)
                txt_tmp = {};
                for j=1:size(mat, 2)
                    txt_tmp{j} = sprintf(spec, scale.*mat(i, j));
                end
                txt{i} = strjoin(txt_tmp, ' , ');
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
        
        function ax = get_figure(name, n)
            % Create a figure and return handle.
            
            fig = findobj( 'Type', 'Figure', 'Name', name);
            if isempty(fig)
                fig = figure('Name', name);
                for i=1:n
                    subplot(1, n, i);
                end
            else
                assert(length(fig)==1, 'invalid figure')
            end
            
            ax = findall(fig, 'type', 'axes');
            assert(length(ax)==n, 'invalid figure')
        end
        
        function get_print(pad, name, varargin)
            % Print a line with padding.
            
            fprintf([pad, name '\n'], varargin{:})
        end
    end
end