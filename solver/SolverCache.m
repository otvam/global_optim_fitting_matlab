classdef SolverCache < handle
    % Class for caching the error function results.
    %
    %    Caching with a specified tolerance for determining unicity.
    %    With or without vectorized call to the error function.
    %    Limit the maxmimum number of elements in the cache.
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
        x_mat_cache % input cached values
        err_mat_cache % output cached values (error)
        wgt_mat_cache % output cached values (weight)
    end
    
    %% public
    methods (Access = public)
        function self = SolverCache(fct_err, cache)
            % Constructor.
            
            % set data
            self.fct_err = fct_err;
            self.use_cache = cache.use_cache;
            self.vec_cache = cache.vec_cache;
            self.n_cache = cache.n_cache;
            self.tol_cache = cache.tol_cache;
            
            % init the cache
            self.get_clear();
        end
        
        function [err_mat, wgt_mat] = get_eval(self, x_mat)
            % Evaluate the error function.
                        
            if self.use_cache==true
                [err_mat, wgt_mat] = self.get_eval_cache(x_mat);
            else
                [err_mat, wgt_mat] = self.get_eval_fct(x_mat);
            end
        end
        
        function get_clear(self)
            % Clear and reset the cache.
            
            self.i_cache = 0;
            self.x_mat_cache = [];
            self.err_mat_cache = [];
            self.wgt_mat_cache = [];
        end
    end
    
    %% private api
    methods (Access = private)
        function [err_mat, wgt_mat] = get_eval_fct(self, x_mat)
            % Call the error function (vectorized or not).
            
            if self.vec_cache
                [err_mat, wgt_mat] = self.fct_err(x_mat);
            else
                for i=1:size(x_mat, 2)
                    [err_mat(:,i), wgt_mat(:,i)] = self.fct_err(x_mat(:,i));
                end
            end
        end
        
        function [err_mat, wgt_mat] = get_eval_cache(self, x_mat)
            % Get the error function through the cache.
            
            % find the element in the cache (NaN if not found)
            for i=1:size(x_mat, 2)
                idx_select(i) = self.get_find(x_mat(:,i));
            end
            
            % get the function values
            [err_mat, wgt_mat] = self.get_value(x_mat, idx_select);
            
            % update the cache
            self.get_update(x_mat, err_mat, wgt_mat, idx_select);
        end
        
        function idx_select = get_find(self, x)
            % Find if an element is in the cache.
            
            if self.i_cache==0
                idx_select = NaN;
            else
                % get the deviation between cache and called value
                err = x-self.x_mat_cache;
                                                
                % get the norm of the error
                err = vecnorm(err, 2, 1);
                
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
        
        function [err_mat, wgt_mat] = get_value(self, x_mat, idx_select)
            % Get the function values (from cache or evaluate).
            
            % call the function if not in the cache
            idx_compute = isnan(idx_select);
            if nnz(idx_compute)==0
                err_mat_compute = [];
                wgt_mat_compute = [];
            else
                [err_mat_compute, wgt_mat_compute] = self.get_eval_fct(x_mat(:,idx_compute));
            end
            
            % get from cache
            idx_cache = isfinite(idx_select);
            idx_access = idx_select(idx_cache);
            err_mat_cache_tmp = self.err_mat_cache(:,idx_access);
            wgt_mat_cache_tmp = self.wgt_mat_cache(:,idx_access);
            
            % assemble the evaluated and cache values
            err_mat(:,idx_compute) = err_mat_compute;
            wgt_mat(:,idx_compute) = wgt_mat_compute;
            err_mat(:,idx_cache) = err_mat_cache_tmp;
            wgt_mat(:,idx_cache) = wgt_mat_cache_tmp;
        end
        
        function get_update(self, x_mat, err_mat, wgt_mat, idx_select)
            % Update the cache.
            
            % remove the used points from the cache
            idx_remove = isfinite(idx_select);
            idx_remove = idx_select(idx_remove);
            self.x_mat_cache(:,idx_remove) = [];
            self.err_mat_cache(:,idx_remove) = [];
            self.wgt_mat_cache(:,idx_remove) = [];
            self.i_cache = self.i_cache-length(idx_remove);
            
            % add the points at the top of the cache
            self.x_mat_cache = [self.x_mat_cache x_mat];
            self.err_mat_cache = [self.err_mat_cache err_mat];
            self.wgt_mat_cache = [self.wgt_mat_cache wgt_mat];
            self.i_cache = self.i_cache+length(idx_select);
            
            % limit cache size
            if self.i_cache>self.n_cache
                idx = 1:(self.i_cache-self.n_cache);
                self.i_cache = self.n_cache;
                
                self.x_mat_cache(:,idx) = [];
                self.err_mat_cache(:,idx) = [];
                self.wgt_mat_cache(:,idx) = [];
            end
        end
    end
end