classdef SolverUtils < handle
    % Static class with various utils.
    %
    %    Computing error metrics.
    %    Scaling of variables.
    %
    %    Thomas Guillod.
    %    2021-2022 - BSD License.

    %% error
    methods(Static, Access = public)
        function err = get_norm(err_mat, wgt_mat, norm)
            % Get the norm on a error vector with specified weights.
            
            % compute the norm
            if isfinite(norm)
                err_wgt_mat = abs(err_mat).*(wgt_mat.^(1./norm));
                n_elem = sum(wgt_mat, 2);
                err = (sum(err_wgt_mat.^norm, 2)./n_elem).^(1./norm);
            elseif isinf(norm)
                err = max(abs(err_mat), [], 2);
            else
                error('invalid norm')
            end
            
            % check for invalid data
            idx = any(isfinite(err_mat)==false, 2);
            err(idx) = NaN;
        end
        
        function err = get_percentile(err_mat, percentile)
            % Get the percentile on a error vector.
            
            % compute the specified percentile
            err = quantile(abs(err_mat), percentile, 2);
            
            % check for invalid data
            idx = any(isfinite(err_mat)==false, 2);
            err(idx) = NaN;
        end
    end
    
    %% scaling
    methods(Static, Access = public)
        function [x_unclamp, lb_unclamp, ub_unclamp] = get_var_unclamp(x, clamp_bnd, lb, ub)
            % Transform bounded a variable into a unconstrained variable with sine transformation.
            
            if clamp_bnd==true
                % the new variable is unconstrained
                lb_unclamp = -Inf;
                ub_unclamp = +Inf;
                
                if isinf(lb)&&isinf(ub)
                    % do nothing as the variable is already unconstrained
                    x_unclamp = x;
                elseif isinf(lb)
                    % quadratic transformation for single-sided bound
                    x_unclamp = sqrt(x-lb);
                elseif isinf(ub)
                    % quadratic transformation for single-sided bound
                    x_unclamp = sqrt(ub-x);
                else
                    % sine transformation for double-sided bounds
                    x_unclamp = 2.*(x-lb)./(ub-lb)-1;
                    x_unclamp = asin(x_unclamp);

                    % shift to avoid numerical issue around zero
                    x_unclamp = 2.*pi+x_unclamp;
                end
            else
                x_unclamp = x;
                lb_unclamp = lb;
                ub_unclamp = ub;
            end
            
            % avoid numerical issue at the bounds
            if x_unclamp<=lb_unclamp
                x_unclamp = lb_unclamp+eps;
            end
            if x_unclamp>=ub_unclamp
                x_unclamp = ub_unclamp-eps;
            end
        end
        
        function x = get_var_clamp(x_unclamp, clamp_bnd, lb, ub)
            % Transform a unconstrained variable into a bounded variable with sine transformation.
            
            if clamp_bnd==true
                if isinf(lb)&&isinf(ub)
                    % do nothing as the variable is already unconstrained
                    x = x_unclamp;
                elseif isinf(lb)
                    % quadratic transformation for single-sided bound
                    x = lb+x_unclamp.^2;
                elseif isinf(ub)
                    % quadratic transformation for single-sided bound
                    x = ub-x_unclamp.^2;
                else
                    % sine transformation for double-sided bounds
                    x = (sin(x_unclamp)+1)./2;
                    x = x.*(ub-lb)+lb;
                end
            else
                x = x_unclamp;
            end
        end

        function [x0_scale, lb_scale, ub_scale] = get_var_scale(x0, lb, ub, scale, norm)
            % Scale a variable (transformation and normalization).
            
            % get the variable transformation
            switch scale
                case 'lin'
                    fct_scale = @(x) x;
                case 'sqrt'
                    fct_scale = @(x) sqrt(x);
                case 'quad'
                    fct_scale = @(x) x.^2;
                case 'log'
                    fct_scale = @(x) log10(x);
                otherwise
                    error('invalid data')
            end
            
            % transform the variable and the bounds
            x0_scale = fct_scale(x0);
            lb_scale = fct_scale(lb);
            ub_scale = fct_scale(ub);
            
            % normalize if required
            if norm==true
                fct = @(x) (x-lb_scale)./(ub_scale-lb_scale);
                
                x0_scale = fct(x0_scale);
                lb_scale = fct(lb_scale);
                ub_scale = fct(ub_scale);
            end
        end
        
        function x = get_var_unscale(x_scale, lb, ub, scale, norm)
            % Unscale a variable (transformation and normalization).
            
            % get the variable transformation
            switch scale
                case 'lin'
                    fct_scale = @(x) x;
                    fct_unscale = @(x) x;
                case 'sqrt'
                    fct_scale = @(x) sqrt(x);
                    fct_unscale = @(x) x.^2;
                case 'quad'
                    fct_scale = @(x) x.^2;
                    fct_unscale = @(x) sqrt(x);
                case 'log'
                    fct_scale = @(x) log10(x);
                    fct_unscale = @(x) 10.^x;
                otherwise
                    error('invalid data')
            end

            % denormalize if required
            if norm==true
                lb_scale = fct_scale(lb);
                ub_scale = fct_scale(ub);
                x_scale = lb_scale+x_scale.*(ub_scale-lb_scale);
            end
            
            % revert the variable transformation
            x = fct_unscale(x_scale);
        end
    end
end