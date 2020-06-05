classdef scaling_factor
    properties
        name;
        value;
    end
    
    methods
        function obj = scaling_factor(varargin)
            if(nargin == 0)
                par = "*";
                val = NaN;
            elseif(nargin == 1)
                v1 = varargin{1};
                if(iscell(v1))
                    v1 = struct(v1{:});
                end
                
                if(isstruct(v1))
                    par = fields(v1);
                    n = length(par);
                    val = zeros(n, 1);
                    for idx = 1:n
                        val(idx) = v1.(par{idx});
                    end
                else
                    error('it should be a structure or a cell.');
                end
            elseif(nargin == 2)
                par = varargin{1};
                val = varargin{2};
            end
            
            if(numel(par) ~= numel(val))
                error("Name and value should have the same number of elements.");
            end
            
            if(isrow(val))
                val = val';
            end
            
            obj.name = par;
            obj.value = val;
        end
        
        function disp(obj)
            if(isnumeric([obj.name{:}]))
                formatSpec = "%d:\t%s\n";
            else
                formatSpec = "%s:\t%.3e\n";
            end
            
            for idx = 1:numel(obj.name)
                fprintf(formatSpec, obj.name{idx}, obj.value(idx));
            end
        end
        
        function val = length(obj)
            val = length(obj.value);
        end
        
%         function val = numel(obj) 
%         % removed because it was messing up subsref
%             val = numel(obj.value);
%         end
        
        function val = size(obj)
            val = [length(obj) 1];
        end
        
        function obj_struct = to_struct(obj)
            obj_cell = obj.to_cell();
            obj_struct = struct(obj_cell{:});
        end
        
        function obj_cell = to_cell(obj)
            obj_cell = cell(0);
            for idx = 1:obj.length()
                obj_cell = [ obj_cell(:)', ...
                            {obj.name{idx}, obj.value(idx)}];
            end
        end
        
        function gamma_sub = ends_with(obj, pattern)
            gm_struct = obj.to_struct();
            field_gm  = fields(gm_struct);
            flag      = endsWith(field_gm, pattern);
            field_gm  = field_gm(flag);
            
            gm = cell(0);
            for idx = 1:length(field_gm)
                new_field = field_gm{idx}(1:end-1);
                if(endsWith(new_field, '_'))
                    new_field = new_field(1:end-1);
                end
                
                gm = [gm(:)', ...
                      {new_field, gm_struct.(field_gm{idx})}];
            end
            gamma_sub = scaling_factor(gm);
        end
        
        function obj_new = update(obj, gamma)
            if(~isempty(gamma))
                obj_struct = obj.to_struct();
                gam_struct = gamma.to_struct();

                fld_gam = fields(gam_struct);
                
                out = obj_struct;
                
                for idx = 1:length(fld_gam)
                    f = fld_gam{idx};
                    if(isfield(obj_struct, f))
                        out.(f) = gam_struct.(f);
                    end
                end
                obj_new = scaling_factor(out);
            else
                obj_new = obj;
            end
                
        end
        
        function obj_new = add_factor(obj, varargin)
            obj_cell = obj.to_cell();
            obj_new = scaling_factor([obj_cell(:)', varargin(:)']);
        end
        
        function some_name = name_contains(obj, str)
            idx = contains(obj.name, str);
            some_name = string(obj.name(idx));
        end
        
    end
    
    %% Set methods:
    methods
        function obj = subsasgn(obj, par, varargin)
            if isequal(obj,[])
                % obj = ClassName.empty;
            end
            
            if(isrow(varargin{:}))
                varargin{:} = varargin{:}';
            end
            
            switch(par(1).type)
                case '.'
                    if length(par) == 1
                        obj.(par.subs) = varargin{:};
                    elseif((length(par) == 2) && (strcmp(par(2).type, '()')))
                        obj.(par(1).subs)(par(2).subs{:}) = varargin{:};
                    else
                        obj = builtin('subsasgn', obj, par, varargin{:});
                    end
                case '()'
                    par_idx = par.subs{:};

                    if(isnumeric(par_idx))
                        if((min(par_idx) >= 1) && (max(par_idx) <= length(obj)))
                            idx = par_idx;
                        else
                            error("Indexes out of range.");
                        end
                    elseif(isstring(par_idx) || ischar(par_idx))
                        par_idx(par_idx == "*") = [];
                        idx = ismember(obj.name, par_idx);
                    end
                    
                    obj.value(idx) = varargin{:}';
                otherwise
                    error("Operator [%s] not implemented yet.", par(1).type);
            end
        end
    end
    
    %% Get methods:
    methods
        function val = get.value(obj)
            % based on: (10/02/2020)
            % https://www.mathworks.com/matlabcentral/answers/337824-why-jsonencode-with-containers-map-doesn-t-preserve-order-of-key-value-pairs
            val = obj.value;
        end
        
        function val = get.name(obj)
            val = obj.name;
        end
        
        function varargout = subsref(obj, par)
            
            switch(par(1).type)
                case '.'
                    if(length(par) == 1)
                        val = obj.(par.subs);
                    elseif((length(par) == 2) && (strcmp(par(2).type,'()')))
                        val = obj.(par(1).subs)(par(2).subs{:});
                    else
                        val = builtin('subsref', obj, par);
                    end
                case '()'
                    par_idx = par.subs{:};

                    if(isnumeric(par_idx))
                        if((min(par_idx) >= 1) && (max(par_idx) <= length(obj)))
                            idx = par_idx;
                        else
                            error("Indexes out of range.");
                        end
                    elseif(islogical(par_idx))
                        if((numel(par_idx) > 1) && (numel(par_idx) <= length(obj)))
                            idx = par_idx;
                        else
                            error("Indexes out of range.");
                        end
                    elseif(isstring(par_idx) || ischar(par_idx))
                        idx = ismember(obj.name, par_idx);
                    end

                    val = obj.value(idx);
                otherwise
                    error("Operator [%s] not implemented yet.", par(1).type);
            end
            
            [varargout{1:nargout}] = val;
        end
        
        function idx = end(obj, k, n)
            szd = size(obj);
            if k < n
                idx = szd(k);
            else
                idx = prod(szd(k:end));
            end
        end
    end
    
    %% Static method:
    methods(Static)
        function gamma = generate_gamma(n, N)
            %GENERATE_GAMMA generates an array of values from power(10, -n) to 1 with N
            % elements.
            %
            
            exp_array = -n:0;
            specific = 10.^exp_array;
            
            p = ceil(N/n);
            gm_lin = zeros(n, p);
            gm_log = zeros(n, p);
            
            for idx = 1:n
                gm_lin(idx, :) = linspace(specific(idx) , specific(idx + 1) , p);
                gm_log(idx, :) = logspace(exp_array(idx), exp_array(idx + 1), p);
            end
            
            gm_lin = reshape(gm_lin, 1, n*p);
            gm_log = reshape(gm_log, 1, n*p);
            
            gm = sort([gm_lin, gm_log]);
            gamma = unique(gm(1:2:end));
        end
    end
    
end
