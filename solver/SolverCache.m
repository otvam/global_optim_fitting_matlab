classdef SolverCache < handle
    % Class for caching the error function results.
    %
    %    Caching with a specified tolerance for determining unicity.
    %    With or without vectorized call to the error function.
    %    Maxmimum number of elements in the cache.
    %
    %    Thomas Guillod.
    %    2021-2022 - BSD License.
    
    %% properties
    properties (SetAccess = private, GetAccess = private)
        fct_err % error function to be cached
        use_cache % use (or not) the cache for the function
        vec_cache % allow (or not) vectorized call to the function
        n_cache % maximum number of elements in the cache
        tol_cache % tolerance for determining unicity
        
        i_cache % number of element currently in the cache
        x_val % input cached values
        err_val % output cached values (error)
        wgt_val % output cached values (weight)
    end
    
    %% public
    methods (Access = public)
        function self = SolverCache(fct_err, use_cache, vec_cache, n_cache, tol_cache)
            % Constructor.
            
            % set data
            self.fct_err = fct_err;
            self.use_cache = use_cache;
            self.vec_cache = vec_cache;
            self.n_cache = n_cache;
            self.tol_cache = tol_cache;
            
            % init the cache
            self.get_clear();
        end
        
        function [err, wgt] = get_eval(self, x)
            % Evaluate the error function.
            
            if self.use_cache==true
                [err, wgt] = self.get_eval_cache(x);
            else
                [err, wgt] = self.get_eval_fct(x);
            end
        end
        
        function get_clear(self)
            % Clear and reset the cache.
            
            self.i_cache = 0;
            self.x_val = [];
            self.err_val = [];
            self.wgt_val = [];
        end
    end
    
    %% private api
    methods (Access = private)
        function [err, wgt] = get_eval_fct(self, x)
            % Call the error function (vectorized or not).
            
            if self.vec_cache
                [err, wgt] = self.fct_err(x);
            else
                for i=1:size(x, 1)
                    [err(i,:), wgt(i,:)] = self.fct_err(x(i,:));
                end
            end
        end
        
        function [err, wgt] = get_eval_cache(self, x)
            % Get the error function through the cache.
            
            % find the element in the cache (NaN if not found)
            for i=1:size(x, 1)
                idx_select(i) = self.get_find(x(i,:));
            end
            
            % get the function values
            [err, wgt] = self.get_value(x, idx_select);
            
            % update the cache
            self.get_update(x, err, wgt, idx_select);
        end
        
        function idx_select = get_find(self, x)
            % Find if an element is in the cache.
            
            if self.i_cache==0
                idx_select = NaN;
            else
                % get the deviation between cache and called value
                err = x-self.x_val;
                
                % NaN is considered as unique
                idx = isnan(x)&isnan(self.x_val);
                err(idx) = 0;
                
                % get the norm of the error
                err = vecnorm(err, 2, 2);
                
                % find the elements matching the tolerance
                idx_remove = find(err<self.tol_cache);
                
                if isempty(idx_remove)
                    idx_select = NaN;
                else
                    % take the best match
                    [~, idx_select] = min(err(idx_remove));
                    
                    % find the global index
                    idx_select = idx_remove(idx_select);
                end
            end
        end
        
        function [err, wgt] = get_value(self, x, idx_select)
            % Get the function values (from cache or evaluate).
            
            % call the function if not in the cache
            idx_compute = isnan(idx_select);
            if nnz(idx_compute)==0
                err_compute = [];
                wgt_compute = [];
            else
                [err_compute, wgt_compute] = self.get_eval_fct(x(idx_compute,:));
            end
            
            % get from cache
            idx_cache = isfinite(idx_select);
            idx_access = idx_select(idx_cache);
            err_cache = self.err_val(idx_access,:);
            wgt_cache = self.wgt_val(idx_access,:);
            
            % assemble the evaluated and cache values
            err(idx_compute,:) = err_compute;
            wgt(idx_compute,:) = wgt_compute;
            err(idx_cache,:) = err_cache;
            wgt(idx_cache,:) = wgt_cache;
        end
        
        function get_update(self, x, err, wgt, idx_select)
            % Update the cache.
            
            % remove the used points from the cache
            idx_remove = isfinite(idx_select);
            idx_remove = idx_select(idx_remove);
            self.x_val(idx_remove,:) = [];
            self.err_val(idx_remove,:) = [];
            self.wgt_val(idx_remove,:) = [];
            self.i_cache = self.i_cache-length(idx_remove);
            
            % add the points at the top of the cache
            self.x_val = [self.x_val ; x];
            self.err_val = [self.err_val ; err];
            self.wgt_val = [self.wgt_val ; wgt];
            self.i_cache = self.i_cache+length(idx_select);
            
            % limit cache size
            if self.i_cache>self.n_cache
                idx = 1:(self.i_cache-self.n_cache);
                self.i_cache = self.n_cache;
                
                self.x_val(idx,:) = [];
                self.err_val(idx,:) = [];
                self.wgt_val(idx,:) = [];
            end
        end
    end
end