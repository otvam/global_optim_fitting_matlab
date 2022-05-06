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
        tol_bound % tolerance for determining if a value is close to a bound
        
        x0_scale % initial values (scaled)
        lb_scale % lower bounds (scaled)
        ub_scale % upper bounds (scaled)
        tol_scale % tolerance on the bounds (scaled)
    end
    
    %% public
    methods (Access = public)
        function self = SolverVar(var_opt, var_fix, tol_bound)
            % Constructor.

            % set data
            self.var_opt = var_opt;
            self.var_fix = var_fix;
            self.tol_bound = tol_bound;
            
            % get the scaled data
            for i=1:length(self.var_opt)
                [self.x0_scale(:,i), self.lb_scale(i), self.ub_scale(i), self.tol_scale(i)] = SolverVar.get_init(self.var_opt{i});
            end
        end
        
        function x0_scale = get_x0_scale(self)
            % Get the initial scaled values.
            
            x0_scale = self.x0_scale;
        end
        
        function [n_pts, param, bnd, is_bound] = get_param(self, x_scale)
            % Extract the parameter structure from a raw matrix.
            %    - unscale the values (transformation and normalization)
            %    - check if the values are closed to the bounds
            %    - assign the results in structs
            
            % handle the fitting variables
            for i=1:length(self.var_opt)
                [name, idx, x_tmp, is_bound_tmp] = SolverVar.get_param_opt(x_scale(:,i), self.var_opt{i}, self.lb_scale(i), self.ub_scale(i), self.tol_bound);
                
                param.(name)(idx,:) = x_tmp;
                bnd.(name)(idx,:) = is_bound_tmp;
                is_bound_opt(i,:) = is_bound_tmp;
            end
            
            % add the fixed variables
            n_rep = size(x_scale, 1);
            for i=1:length(self.var_fix)
                [name, idx, x_tmp, is_bound_tmp] = SolverVar.get_param_fix(self.var_fix{i}, n_rep);
                
                param.(name)(idx,:) = x_tmp;
                bnd.(name)(idx,:) = is_bound_tmp;
                is_bound_fix(i,:) = is_bound_tmp;
            end
                        
            % order (otherwise random)
            param = orderfields(param);
            bnd = orderfields(bnd);

            % get the number of parameter combinations
            n_pts = size(x_scale, 1);

            % check if any parameters are close to the bounds
            is_bound = all(is_bound_opt, 1)&all(is_bound_fix, 1);
        end
        
        function [x_unclamp, lb_unclamp, ub_unclamp] = get_unclamp(self, x_scale, clamp_bnd)
            % Transform bounded variables into unconstrained variables with sine transformation.
            
            for i=1:size(x_scale, 2)
                [x_unclamp(:,i), lb_unclamp(i), ub_unclamp(i)] = SolverUtils.get_var_unclamp(x_scale(:,i), clamp_bnd, self.lb_scale(i), self.ub_scale(i));
            end
        end
        
        function x_scale = get_clamp(self, x_unclamp, clamp_bnd)
            % Transform unconstrained variables into bounded variables with sine transformation.
            
            for i=1:size(x_unclamp, 2)
                x_scale(:,i) = SolverUtils.get_var_clamp(x_unclamp(:,i), clamp_bnd, self.lb_scale(i), self.ub_scale(i));
            end
        end
    end
    
    %% private static api
    methods(Static, Access = private)
        function [x0_scale, lb_scale, ub_scale, tol_scale] = get_init(var)
            % Scale a variable (bounds, transformation, and normalization).
            
            % extract
            x0 = var.x0;
            lb = var.lb;
            ub = var.ub;
            tol_bnd = var.tol_bnd;
            trf = var.trf;
            norm = var.norm;
            
            % scale the variable
            x0_scale =  SolverUtils.get_var_scale(x0, lb, ub, trf, norm);
            lb_scale =  SolverUtils.get_var_scale(lb, lb, ub, trf, norm);
            ub_scale =  SolverUtils.get_var_scale(ub, lb, ub, trf, norm);            
            tol_scale = tol_bnd.*(ub_scale-lb_scale);
        end
        
        function [name, idx, x, is_bound] = get_param_opt(x_scale, var, lb_scale, ub_scale, tol_bound)
            % Unscale a variable (bounds, transformation, and normalization).
            
            % extract
            idx = var.idx;
            name = var.name;
            lb = var.lb;
            ub = var.ub;
            trf = var.trf;
            norm = var.norm;
                        
            % check if any parameters are close to the bounds 
            tol = tol_bound.*(ub_scale-lb_scale);
            lb_scale_tol = lb_scale+tol;
            ub_scale_tol = ub_scale-tol;
            is_bound = (x_scale>lb_scale_tol)&(x_scale<ub_scale_tol);
            
            % unscale the variable
            x = SolverUtils.get_var_unscale(x_scale, lb, ub, trf, norm);
        end
        
        function [name, idx, x, is_bound] = get_param_fix(var, n_rep)
            % Get the values of a fixed parameter.
            
            idx = var.idx;
            name = var.name;
            x = var.x;
            
            x = repmat(x, n_rep, 1);
            is_bound = true(n_rep, 1);
        end
    end
end