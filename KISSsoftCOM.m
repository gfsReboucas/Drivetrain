classdef KISSsoftCOM
    %KISSSOFTCOM Wrapper class to deal with KISSsoft's COM interface. Based
    % on code provided by KISSsoft support team. More details about
    % KISSsoft can be found on [1].
    % [1] https://www.kisssoft.com/en
    %
    
    properties(Access = private)
        COM;
    end
    
    properties
        module;
    end
    
    methods
        function obj = KISSsoftCOM(mod)
            [flag, err] = KISSsoftCOM.is_installed();
            if(~flag)
                error(err.identifier, "%s", err.message);
            end
            
            obj.module = mod;

            obj.COM.SetSilentMode(true);
            obj.COM.GetModule(obj.module, false);
        end
        
        function set_var(obj, name, val)
            if(~obj.check_variable(name))
                error('Variable [%s] does not exist.', upper(name));
            end
            
            flag = obj.COM.SetVar(name, num2str(val, 10));
            if(flag == false)
                error("Could not set variable [%s] to %.3f.", name, val);
            end
        end
        
        function val = get_var(obj, name)
            if(~obj.is_module_loaded())
                warning('No module loaded');
            end
            
            val = obj.COM.GetVar(name);
            if(isnan(val))
                error('Variable [%s] does not exist or cannot be accessed.', upper(name));
            end
            val = str2double(val);
        end
        
        function val = get_version(obj)
            val = obj.COM.GetKsoftVersion();
        end
        
        function val = calculate(obj)
            val = true;
            if(~obj.is_module_loaded())
                warning('No module loaded');
                val = false;
            end
            
            if(~obj.COM.CalculateRetVal())
                error('Error(s) during calculation');
            else
                val = true;
            end
        end
                
        function load_file(obj, file_name)
            obj.COM.LoadFile(file_name);
        end
        
        function val = check_variable(obj, name)
            % check if variable exist.
            val = true;
            
            if(~obj.is_module_loaded())
                warning('No module loaded');
                val = false;
            end
            
            if(isnan(obj.COM.GetVar(name)))
                val = false;
            else
                val = true;
            end
        end
        
        function val = is_module_loaded(obj)
            % checks if a module is loaded.
            val = obj.COM.IsActive();
        end
        
        function save_file(obj, file_name)
            if(~obj.is_module_loaded())
                warning('No module loaded');
            end
            
            obj.COM.SaveFile(file_name);
        end
        
        function write_report(obj, template, name, show_flag, format)
            obj.COM.ReportWithParameters(template, name, show_flag, format);
        end
        
        function report(obj, flag)
            obj.COM.Report(flag);
        end
        
        function delete(obj)
            obj.COM.ReleaseModule();
        end
    end
    
    methods(Static)
        function [val, err] = is_installed()
            val = true;
            err = [];
            try
                actxserver('KISSsoftCOM.KISSsoft');
            catch err
                val = false;
            end
        end
    end
end
