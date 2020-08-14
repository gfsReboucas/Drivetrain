classdef KISSsoftCOM
    %KISSSOFTCOM Wrapper class to deal with KISSsoft's COM interface. Based
    % on code provided by KISSsoft support team. More details about
    % KISSsoft can be found on [1]. Run KISSsoftCOM.view_COM_methods() to
    % view all the methods provided by KISSsoft's COM interface. More info
    % about it at [2].
    %
    % [1] https://www.kisssoft.com/en
    % [2] https://mathworks.com/help/matlab/ref/methodsview.html
    %
    
    properties(Access = private)
        COM;
    end
    
    properties
        module;
    end
    
    methods
        function obj = KISSsoftCOM(mod)
            [flag, msg] = KISSsoftCOM.is_installed();
            if(flag)
                obj.COM = actxserver('KISSsoftCOM.KISSsoft');
            else
                error(msg.identifier, "%s", msg.message);
            end
            
            obj.module = mod;

            obj.COM.SetSilentMode(false);
            obj.COM.GetModule(obj.module, false);
        end
        
        function show_UI(obj)
            obj.COM.SetSilentMode(false);
            obj.COM.ShowInterface(true);
        end
        
        function set_var(obj, name, val)
            if(~obj.check_variable(name))
                error('Variable [%s] does not exist.', upper(name));
            end
            
            obj.COM.SetVar(name, num2str(val, 10));
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
        
        function calculate(obj)
            obj.COM.Calculate();
        end
        
        function val = check_calculations(obj)
            flag_mod = obj.is_module_loaded();
            flag_calc = obj.COM.CalculateRetVal();
            
            if(flag_mod == false)
                val = false;
                warning('No module loaded');
            else
                val = true;
            end
            
            
            if(flag_calc == false)
                val = false;
                warning('Error(s) during calculation');
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
            
            try
                obj.COM.GetVar(name);
            catch err
                val = false;
                warning(err.identifier, '%sVariable [%s] undefined.', err.message, upper(name));
            end
            
        end
        
        function val = is_module_loaded(obj)
            % checks if a module is loaded.
            val = obj.COM.isActive();
        end
        
        function save_file(obj, file_name)
            if(~obj.is_module_loaded())
                warning('No module loaded');
            end
            
            obj.COM.SaveFile([file_name, '.',strrep(obj.module, '0', '')]);
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
        function [val, msg] = is_installed()
            val = true;
            msg = [];
            try
                actxserver('KISSsoftCOM.KISSsoft');
            catch msg
                val = false;
            end
        end
        
        function view_COM_methods()
            [flag, msg] = KISSsoftCOM.is_installed();
            if(flag)
                methodsview(actxserver('KISSsoftCOM.KISSsoft'));
            else
                warning(msg.identifier, "%s", msg.message);
            end
        end
        
        function example_02()
            %EXAMPLE_02 from KUM
            %
            
            ex02 = KISSsoftCOM('Z012');
            version = ex02.get_version();
            version = strrep(version, '/', '-');            
            
            std_file = 'CylGearPair 1 (spur gear).Z12';
            file_name = sprintf('C:\\Program Files (x86)\\KISSsoft %s\\example\\%s', version, std_file);
            
            ex02.load_file(file_name);
            
            ex02.check_calculations();
            ex02.calculate();
            val = ex02.get_var('ZPP[0].Fuss.SFnorm');
            disp(val);
        end
        
        function example_03()
            %EXAMPLE_03 from KUM
            %
            
            ex03 = KISSsoftCOM('Z012');
            version = ex03.get_version();
            version = strrep(version, '/', '-');            
            
            std_file = 'CylGearPair 1 (spur gear).Z12';
            file_name = sprintf('C:\\Program Files (x86)\\KISSsoft %s\\example\\%s', version, std_file);
            
            ex03.load_file(file_name);
            
            ex03.check_calculations();
            ex03.calculate();
            
%             fprintf('ZF_original = %.3f\n', ex03.get_var("c"));
            
            fprintf('Current profile shift gear 1 : %.3f\n', ex03.get_var("ZR[0].x.nul"));
            fprintf('Current tooth root safety gear 1: %.3f\n', ex03.get_var("ZPP[0].Fuss.SFnorm"));
            
            ex03.set_var('ZR[0].x.nul', 0.1);
            
            ex03.check_calculations();
            ex03.calculate();
            
            fprintf('New profile shift gear 1 : %.3f\n', ex03.get_var("ZR[0].x.nul"));
            fprintf('New tooth root safety gear 1: %.3f\n', ex03.get_var("ZPP[0].Fuss.SFnorm"));

        end
        
        
        
    end
end