function out = process_varargin(standard, updated)
    standard = struct(standard{:});
    
    if(~exist('updated', "var"))
        updated = cell(0);
    end
    
    if(~isempty(updated))
        arg_in = struct(updated{:});
        out = standard;
        
        field_v = fields(arg_in);
        
        lv = length(field_v);
        
        for idx = 1:lv
            f = field_v{idx};
            if(isfield(standard, f))
                out.(f) = arg_in.(f);
            end
        end
    else
        out = standard;
    end
    
end
