classdef SolverVar < handle
    % Class for managing the variables.
    %
    %    Abstraction layer (high level parameter structure vs. raw matrix).
    %    Variable transformation (linear, quadratic, logarithmic).
    %    Variable normalization (between zero and one).
    %    Transform bounded variables to uncontrained variable with sine transformation.
    %
    %    Thomas Guillod.
    %    2021-2022 - BSD License.
    
    %% properties
    properties (SetAccess = private, GetAccess = private)
        var_opt % description of the parameters to be fitted
        var_fix % description of the parameters with fixed values
    end
    
    %% public
    methods (Access = public)
        function self = SolverVar(var_opt, var_fix)
            % Constructor.

            % set data
            self.var_opt = var_opt;
            self.var_fix = var_fix;
        end
        
        function [n_pts, param] = get_init(self)
            % Get the initial values.
            
            % handle the fitting variables
            for i=1:length(self.var_opt)
                x0 = self.var_opt{i}.x0;
                idx = self.var_opt{i}.idx;
                name = self.var_opt{i}.name;
                
                n_pts(i) = length(x0);
                param.(name)(idx,:) = x0;
            end
            
            % get number of points
            n_pts = unique(n_pts);
            assert(length(n_pts)==1, 'invalid data')
            
            % handle the fixed variables
            for i=1:length(self.var_fix)
                x0 = self.var_fix{i}.x0;
                idx = self.var_fix{i}.idx;
                name = self.var_fix{i}.name;
                
                param.(name)(idx,:) = repmat(x0, n_pts, 1);
            end
        end
        
        function [n_var, x_scale, lb_scale, ub_scale] = get_scale(self, n_pts, param, clamp_bnd)
            % Extract the raw matrix from the parameter structure.
            %    - scale the values (transformation, normalization, and clamping)
            %    - handle the bounds
            %    - assign the results in matrices

            % handle the fitting variables
            for i=1:length(self.var_opt)
                % extract
                idx = self.var_opt{i}.idx;
                name = self.var_opt{i}.name;
                lb = self.var_opt{i}.lb;
                ub = self.var_opt{i}.ub;
                trf = self.var_opt{i}.trf;
                norm = self.var_opt{i}.norm;
                
                % get the bounds
                lb_trf = SolverUtils.get_var_trf(lb, trf, false);
                ub_trf = SolverUtils.get_var_trf(ub, trf, false);
                lb_norm = SolverUtils.get_var_norm(lb_trf, lb_trf, ub_trf, norm, false);
                ub_norm = SolverUtils.get_var_norm(ub_trf, lb_trf, ub_trf, norm, false);

                % get vector
                x_tmp = param.(name)(idx,:);
                assert(length(x_tmp)==n_pts, 'invalid data')
                
                % variable scaling
                x_trf_tmp = SolverUtils.get_var_trf(x_tmp, trf, false);
                x_norm_tmp = SolverUtils.get_var_norm(x_trf_tmp, lb_trf, ub_trf, norm, false);
                [x_umclamp_tmp, lb_unclamp_tmp, ub_unclamp_tmp] = SolverUtils.get_var_unclamp(x_norm_tmp, lb_norm, ub_norm, clamp_bnd);
                
                % assign
                x_scale(i,:) = x_umclamp_tmp;
                lb_scale(i,:) = lb_unclamp_tmp;
                ub_scale(i,:) = ub_unclamp_tmp;
            end
                                    
            % check the data
            n_var = length(self.var_opt);
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
            is_bound = true;
            for i=1:length(self.var_opt)
                % extract
                idx = self.var_opt{i}.idx;
                name = self.var_opt{i}.name;
                lb = self.var_opt{i}.lb;
                ub = self.var_opt{i}.ub;
                trf = self.var_opt{i}.trf;
                norm = self.var_opt{i}.norm;
                tol_bnd = self.var_opt{i}.tol_bnd;
                
                                % get the bounds
                lb_trf = SolverUtils.get_var_trf(lb, trf, false);
                ub_trf = SolverUtils.get_var_trf(ub, trf, false);
                lb_norm = SolverUtils.get_var_norm(lb_trf, lb_trf, ub_trf, norm, false);
                ub_norm = SolverUtils.get_var_norm(ub_trf, lb_trf, ub_trf, norm, false);
                
                                % get vector
                x_scale_tmp = x_scale(i,:);
                assert(length(x_scale_tmp)==n_pts, 'invalid data')

                % variable scaling
                x_norm_tmp = SolverUtils.get_var_clamp(x_scale_tmp, lb_norm, ub_norm, clamp_bnd);
                x_trf_tmp = SolverUtils.get_var_norm(x_norm_tmp, lb_trf, ub_trf, norm, true);
                x_tmp = SolverUtils.get_var_trf(x_trf_tmp, trf, true);
                
                % check if any parameters are close to the bounds
                tol = tol_bnd.*(ub_norm-lb_norm);
                lb_norm_tol = lb_norm+tol;
                ub_norm_tol = ub_norm-tol;
                is_bound_tmp = is_bound&(x_norm_tmp>lb_norm_tol)&(x_norm_tmp<ub_norm_tol);
                is_bound = is_bound&is_bound_tmp;           
                
                % assign
                param.(name)(idx,:) = x_tmp;
                bnd.(name)(idx,:) = is_bound_tmp;
            end
                        
            % handle the fixed variables
            for i=1:length(self.var_fix)
                x0 = self.var_fix{i}.x0;
                idx = self.var_fix{i}.idx;
                name = self.var_fix{i}.name;
                
                param.(name)(idx,:) = repmat(x0, n_pts, 1);
                bnd.(name)(idx,:) = true(n_pts, 1);
            end
        end
    end
end