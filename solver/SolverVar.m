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
            
            % init the struct
            param = struct();
            
            % handle the fitting variables
            for i=1:length(self.var_opt)
                [n_pts(i), param] = SolverVar.get_var_opt_init(self.var_opt{i}, param);
            end
            
            % get number of points
            n_pts = unique(n_pts);
            assert(length(n_pts)==1, 'invalid data')
            
            % handle the fixed variables
            for i=1:length(self.var_fix)
                param = SolverVar.get_var_fix(self.var_fix{i}, param);
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
                [x_umclamp_tmp, lb_unclamp_tmp, ub_unclamp_tmp] = SolverVar.get_var_opt_scale(self.var_opt{i}, n_pts, param, clamp_bnd);
                                       
                % assign
                x_scale = [x_scale ; x_umclamp_tmp];
                lb_scale = [lb_scale ; lb_unclamp_tmp];
                ub_scale = [ub_scale ; ub_unclamp_tmp];
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
            
            % init
            param = struct();
            bnd = struct();
            is_bound = true;

            % handle the fitting variables
            pos = 0;
            for i=1:length(self.var_opt)
                [param, bnd, is_bound_tmp, pos] = SolverVar.get_var_opt_unscale(self.var_opt{i}, x_scale, clamp_bnd, param, bnd, pos);
                                
                % check if the bounds are globally respected
                is_bound = is_bound&all(is_bound_tmp, 1);
            end
            assert(size(x_scale, 1)==pos, 'invalid data')
            
            % handle the fixed variables
            for i=1:length(self.var_fix)
                param = SolverVar.get_var_fix(self.var_fix{i}, param);
            end
        end
        
        function err_best = get_scale_err(self, err_mat, wgt_mat)
            % Compute the error metric from the output of the error function.
            %    - from the error matrix and the associated weights
            %    - using the specific method
            
            % extract
            type = self.var_err.type;
            arg = self.var_err.arg;
            
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
    
    %% static private api
    methods (Static, Access = private)
        function [n_pts, param] = get_var_opt_init(var_opt, param)
            % Add an initial value to the parameter struct.
            
            % extract
            x0 = var_opt.x0;
            idx = var_opt.idx;
            name = var_opt.name;
            
            % handle the fitting variables
            n_pts = size(x0, 2);
            if isempty(idx)
                param.(name) = x0;
            else
                assert(size(x0, 1)==length(idx), 'invalid data')
                param.(name)(idx,:) = x0;
            end
        end
        
        function [x_umclamp, lb_unclamp, ub_unclamp] = get_var_opt_scale(var_opt, n_pts, param, clamp_bnd)
            % Scale a fitting variable.
            
            % extract
            name = var_opt.name;
            x0 = var_opt.x0;
            idx = var_opt.idx;
            lb = var_opt.lb;
            ub = var_opt.ub;
            trf = var_opt.trf;
            norm = var_opt.norm;
            
            % get variable size
            n_size = size(x0, 1);
            
            % get the bounds
            lb_trf = SolverUtils.get_var_trf(lb, trf, false);
            ub_trf = SolverUtils.get_var_trf(ub, trf, false);
            lb_norm = SolverUtils.get_var_norm(lb_trf, lb_trf, ub_trf, norm, false);
            ub_norm = SolverUtils.get_var_norm(ub_trf, lb_trf, ub_trf, norm, false);
            
            % get vector
            if isempty(idx)
                x = param.(name);
            else
                assert(length(idx)==n_size, 'invalid data')
                x = param.(name)(idx,:);
            end
            assert(size(x, 2)==n_pts, 'invalid data')
            assert(size(x, 1)==n_size, 'invalid data')
            
            % variable scaling
            x_trf = SolverUtils.get_var_trf(x, trf, false);
            x_norm = SolverUtils.get_var_norm(x_trf, lb_trf, ub_trf, norm, false);
            [x_umclamp, lb_unclamp, ub_unclamp] = SolverUtils.get_var_unclamp(x_norm, lb_norm, ub_norm, clamp_bnd);
            
            % make the bounds the same size as the variable
            lb_unclamp = repmat(lb_unclamp, n_size, 1);
            ub_unclamp = repmat(ub_unclamp, n_size, 1);
        end
        
        function [param, bnd, is_bound, pos] = get_var_opt_unscale(var_opt, x_scale, clamp_bnd, param, bnd, pos)
            % Unscale a fitting variable.
            
            name = var_opt.name;
            x0 = var_opt.x0;
            idx = var_opt.idx;
            lb = var_opt.lb;
            ub = var_opt.ub;
            trf = var_opt.trf;
            norm = var_opt.norm;
            tol_bnd = var_opt.tol_bnd;
            
            % check size
            assert(length(size(x0))==2, 'invalid data')

            % get variable size
            n_size = size(x0, 1);
            n_pts = size(x_scale, 2);

            % get the bounds
            lb_trf = SolverUtils.get_var_trf(lb, trf, false);
            ub_trf = SolverUtils.get_var_trf(ub, trf, false);
            lb_norm = SolverUtils.get_var_norm(lb_trf, lb_trf, ub_trf, norm, false);
            ub_norm = SolverUtils.get_var_norm(ub_trf, lb_trf, ub_trf, norm, false);
            
            % get vector
            idx_vec = (pos+1):(pos+n_size);
            x_scale_tmp = x_scale(idx_vec,:);
            assert(size(x_scale_tmp, 1)==n_size, 'invalid data')
            assert(size(x_scale_tmp, 2)==n_pts, 'invalid data')
            
            % update the index
            pos = pos+n_size;
            
            % variable scaling
            x_norm_tmp = SolverUtils.get_var_clamp(x_scale_tmp, lb_norm, ub_norm, clamp_bnd);
            x_trf_tmp = SolverUtils.get_var_norm(x_norm_tmp, lb_trf, ub_trf, norm, true);
            x_tmp = SolverUtils.get_var_trf(x_trf_tmp, trf, true);
            
            % check if any parameters are close to the bounds
            tol = tol_bnd.*(ub_norm-lb_norm);
            lb_norm_tol = lb_norm+tol;
            ub_norm_tol = ub_norm-tol;
            if isfinite(tol)
                is_bound = (x_norm_tmp>lb_norm_tol)&(x_norm_tmp<ub_norm_tol);
            else
                is_bound = (x_norm_tmp>lb_norm)&(x_norm_tmp<ub_norm);
            end
            
            % assign
            if isempty(idx)
                param.(name) = x_tmp;
                bnd.(name) = is_bound;
            else
                assert(length(idx)==n_size, 'invalid data')
                param.(name)(idx,:) = x_tmp;
                bnd.(name)(idx,:) = is_bound;
            end
        end
        
        function param = get_var_fix(var_fix, param)
            % Add a fixed variable to the parameter struct.
            
            % extract
            x0 = var_fix.x0;
            name = var_fix.name;
            
            % check size
            assert(length(size(x0))==2, 'invalid data')
            
            % handle the fixed variables
            param.(name) = x0;
        end
    end
end