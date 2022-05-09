classdef SolverVar < handle
    % Class for managing the variables.
    %
    %    Abstraction layer (high level parameter structure vs. raw matrix).
    %    Variable transformation (linear, quadratic, logarithmic).
    %    Variable normalization (between zero and one).
    %    Transform bounded variables to uncontrained variable with sine transformation.
    %    Compute the desired error metric from the output of the error function.
    %
    %    Thomas Guillod.
    %    2021-2022 - BSD License.
    
    %% properties
    properties (SetAccess = private, GetAccess = private)
        var_opt % description of the parameters to be fitted
        var_fix % description of the parameters with fixed values
        var_err % description of the used error metric
    end
    
    %% public
    methods (Access = public)
        function self = SolverVar(var_opt, var_fix, var_err)
            % Constructor.
            
            % set data
            self.var_opt = var_opt;
            self.var_fix = var_fix;
            self.var_err = var_err;
        end
        
        function [n_pts, param] = get_init(self)
            % Get the variable initial values.
            
            % handle the fitting variables
            for i=1:length(self.var_opt)
                x0 = self.var_opt{i}.x0;
                name = self.var_opt{i}.name;
                
                n_pts(i) = size(x0, 2);
                param.(name) = x0;
            end
            
            % get number of points
            n_pts = unique(n_pts);
            assert(length(n_pts)==1, 'invalid data')
            
            % handle the fixed variables
            for i=1:length(self.var_fix)
                x0 = self.var_fix{i}.x0;
                name = self.var_fix{i}.name;
                
                assert(size(x0, 2)==1, 'invalid data')
                param.(name) = repmat(x0, 1, n_pts);
            end
        end
        
        function [x_scale, lb_scale, ub_scale] = get_scale(self, n_pts, param, clamp_bnd)
            % Extract the raw matrix from the parameter structure.
            %    - scale the values (transformation, normalization, and clamping)
            %    - handle the bounds
            %    - assign the results in matrices
            
            % initialize the matrices
            x_scale = [];
            lb_scale = [];
            ub_scale = [];

            % handle the fitting variables
            for i=1:length(self.var_opt)
                % extract
                name = self.var_opt{i}.name;
                                x0 = self.var_opt{i}.x0;
                lb = self.var_opt{i}.lb;
                ub = self.var_opt{i}.ub;
                trf = self.var_opt{i}.trf;
                norm = self.var_opt{i}.norm;
                
                % get variable size
                n_size = size(x0, 1);

                % get the bounds
                lb_trf = SolverUtils.get_var_trf(lb, trf, false);
                ub_trf = SolverUtils.get_var_trf(ub, trf, false);
                lb_norm = SolverUtils.get_var_norm(lb_trf, lb_trf, ub_trf, norm, false);
                ub_norm = SolverUtils.get_var_norm(ub_trf, lb_trf, ub_trf, norm, false);
                                
                % get vector
                x_tmp = param.(name);
                assert(size(x_tmp, 2)==n_pts, 'invalid data')
                assert(size(x_tmp, 1)==n_size, 'invalid data')
                
                % variable scaling
                x_trf_tmp = SolverUtils.get_var_trf(x_tmp, trf, false);
                x_norm_tmp = SolverUtils.get_var_norm(x_trf_tmp, lb_trf, ub_trf, norm, false);
                [x_umclamp_tmp, lb_unclamp_tmp, ub_unclamp_tmp] = SolverUtils.get_var_unclamp(x_norm_tmp, lb_norm, ub_norm, clamp_bnd);
                                
                % assign
                x_scale = [x_scale ; x_umclamp_tmp];
                lb_scale = [lb_scale ; repmat(lb_unclamp_tmp, n_size, 1)];
                ub_scale = [ub_scale ; repmat(ub_unclamp_tmp, n_size, 1)];
            end
                        
            % check the data
            assert(size(x_scale, 1)==size(lb_scale, 1), 'invalid data')
            assert(size(x_scale, 1)==size(ub_scale, 1), 'invalid data')
            assert(size(x_scale, 2)==n_pts, 'invalid data')
            assert(size(lb_scale, 2)==1, 'invalid data')
            assert(size(ub_scale, 2)==1, 'invalid data')
        end
        
        function [n_pts, param, bnd, is_bound] = get_unscale(self, x_scale, clamp_bnd)
            % Extract the parameter structure from a raw matrix.
            %    - unscale the values (transformation, normalization, and clamping)
            %    - check if the values are closed to the bounds
            %    - assign the results in structure
            
            % get number of points
            n_pts = size(x_scale, 2);
            
            % handle the fitting variables
            idx = 0;
            is_bound = true;
            for i=1:length(self.var_opt)
                % extract
                name = self.var_opt{i}.name;
                                x0 = self.var_opt{i}.x0;
                lb = self.var_opt{i}.lb;
                ub = self.var_opt{i}.ub;
                trf = self.var_opt{i}.trf;
                norm = self.var_opt{i}.norm;
                tol_bnd = self.var_opt{i}.tol_bnd;
                
                % get variable size
                n_size = size(x0, 1);
                
                % get the bounds
                lb_trf = SolverUtils.get_var_trf(lb, trf, false);
                ub_trf = SolverUtils.get_var_trf(ub, trf, false);
                lb_norm = SolverUtils.get_var_norm(lb_trf, lb_trf, ub_trf, norm, false);
                ub_norm = SolverUtils.get_var_norm(ub_trf, lb_trf, ub_trf, norm, false);
                                
                % get vector
                idx_vec = (idx+1):(idx+n_size);
                x_scale_tmp = x_scale(idx_vec,:);
                assert(size(x_scale_tmp, 1)==n_size, 'invalid data')
                assert(size(x_scale_tmp, 2)==n_pts, 'invalid data')
                
                % update the index
                idx = idx+n_size;
                
                % variable scaling
                x_norm_tmp = SolverUtils.get_var_clamp(x_scale_tmp, lb_norm, ub_norm, clamp_bnd);
                x_trf_tmp = SolverUtils.get_var_norm(x_norm_tmp, lb_trf, ub_trf, norm, true);
                x_tmp = SolverUtils.get_var_trf(x_trf_tmp, trf, true);
                
                % check if any parameters are close to the bounds
                tol = tol_bnd.*(ub_norm-lb_norm);
                lb_norm_tol = lb_norm+tol;
                ub_norm_tol = ub_norm-tol;
                is_bound_tmp = is_bound&(x_norm_tmp>lb_norm_tol)&(x_norm_tmp<ub_norm_tol);
                is_bound = is_bound&all(is_bound_tmp, 1);
                                
                % assign
                param.(name) = x_tmp;
                bnd.(name) = is_bound_tmp;
            end
            assert(size(x_scale, 1)==idx, 'invalid data')
            
            % handle the fixed variables
            for i=1:length(self.var_fix)
                x0 = self.var_fix{i}.x0;
                name = self.var_fix{i}.name;
                
                assert(size(x0, 2)==1, 'invalid data')
                x0_tmp = repmat(x0, 1, n_pts);
                is_bound_tmp = true(size(x0_tmp));

                param.(name) = x0_tmp;
                bnd.(name) = is_bound_tmp;
            end
        end
        
        function [err_best, n_set] = get_scale_err(self, err_mat, wgt_mat)
            % Compute the error metric from the output of the error function.
            %    - from the error matrix and the associated weights
            %    - using the specific method

            % extract
            type = self.var_err.type;
            arg = self.var_err.arg;
            
            % check
            assert(all(size(err_mat)==size(wgt_mat)), 'invalid data')
            n_set = (size(err_mat, 1)+size(wgt_mat, 1))./2;
            
            % get the specified error metric
            switch type
                case 'error'
                    err_best = SolverUtils.get_error(err_mat, wgt_mat, arg);
                case 'norm'
                    err_best = SolverUtils.get_norm(err_mat, wgt_mat, arg);
                case 'percentile'
                    err_best = SolverUtils.get_percentile(err_mat, wgt_mat, arg);
                otherwise
                    error('invalid data')
            end
        end
    end
end