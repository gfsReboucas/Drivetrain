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
        folder;
    end
    
    properties
        file_name;
        module;
    end
    
    methods
        function obj = KISSsoftCOM(name_file)
            [flag, msg] = KISSsoftCOM.is_installed();
            if(flag)
                obj.COM = actxserver('KISSsoftCOM.KISSsoft');
            else
                error(msg.identifier, '%s', msg.message);
            end
            
            obj.load_file(name_file);
            
            obj.COM.SetSilentMode(false);
        end
        
        function val = calculate(obj)
            val = true;
            if(~obj.is_module_loaded())
                warning('No module loaded');
                val = false;
            end
            
            if(~obj.COM.CalculateRetVal())
                warning('Error(s) in the calculations');
                val = false;
            end
            
        end
        
        function val = is_module_loaded(obj)
            % checks if a module is loaded.
            val = obj.COM.isActive();
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
%             obj.check_calculations();
        end
        
        function val = get_var(obj, name)
            if(~obj.is_module_loaded())
                warning('No module loaded');
            end
            
%             obj.check_calculations();
            val = obj.COM.GetVar(name);
            if(isnan(val))
                error('Variable [%s] does not exist or cannot be accessed.', upper(name));
            end
            val = str2double(val);
        end
        
        function val = check_calculations(obj)
            flag_mod = obj.is_module_loaded();
            flag_calc = obj.COM.CalculateRetVal();
            KS_msg = obj.message();
            
            if(flag_mod == false)
                val = false;
                warning('KISSsoftCOM:no_module', 'No module loaded.\n\tKISSsofts''s message:\n%s', KS_msg);
            else
                val = true;
            end
            
            
            if(flag_calc == false)
                val = false;
                warning('KISSsoftCOM:wrong_calc', 'Problem(s) during calculation.\n\tKISSsofts''s message:\n%s', KS_msg);
            else
                val = true;
            end
        end
                
        function val = load_file(obj, name_file)
            
            if(~obj.is_module_loaded())
                obj.set_module(name_file);
            end
            
            % Compare messages before and after loading file
            msg_bef = obj.message();
            obj.COM.LoadFile(name_file);
            msg_aft = obj.message();
            
            if(~strcmpi(msg_bef, msg_aft))
                error('Could not load file [%s].\n\tKISSsoft''s message:\n%s', upper(name_file), msg_aft);
            end
            val = true;
        end
        
        function val = load_example_file(obj, file_name)
            val = obj.load_file([obj.folder, '\example\', file_name]);
        end
        
        function val = check_variable(obj, name)
            % check if variable exist.
            val = true;
            
            if(~obj.is_module_loaded())
                warning('No module loaded');
                val = false;
            end
            
            if(isempty(obj.COM.GetVar(name)))
                warning('KISSsoftCOM:check_var', ...
                        'Variable [%s] undefined.', upper(name));
                val = false;
            end
            
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
        
        function [fold, fname, mod] = set_module(obj, name_file)
%             if(~isfile(name_file))
%                 error('File [%s] does not exist', upper(name_file));
%             end
%             
            % Get folder, file name, and module from file:
            [fold, ...
             fname, ...
             mod] = fileparts(name_file);
            number = mod(3:end);
            if(isnan(str2double(number)))
                error('KISSsoftCOM:set_module', 'Wrong file extension [%s]', upper(mod));
            end
            mod = sprintf('%s0%s', mod(2), number);
            obj.COM.GetModule(mod, false);
        end
        
        function msg = message(obj)
            msg = obj.COM.Message();
            msg = cell2mat(msg');
            if(isempty(msg))
                msg = '';
            end
        end
        
    end
    
    methods(Static)
        function val = get_version()
            val = nan;
            if(KISSsoftCOM.is_installed())
                val = actxserver('KISSsoftCOM.KISSsoft').GetKsoftVersion();
            end
        end
        
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
                warning(msg.identifier, '%s', msg.message);
            end
        end
        
        function example_02()
            %EXAMPLE_02 from KISSsoft's Users Meeting - KUM.
            %
            
            ex02 = KISSsoftCOM('Z012');
            
            std_file = 'CylGearPair 1 (spur gear).Z12';
%             file_name = sprintf('C:\\Program Files\\KISSsoft %s\\example\\%s', version, std_file);
            
            ex02.load_example_file(std_file);
            
            ex02.check_calculations();
            ex02.calculate();
            val = ex02.get_var('ZPP[0].Fuss.SFnorm');
            disp(val);
        end
        
        function example_03()
            %EXAMPLE_03 from KUM
            %
            
            ex03 = KISSsoftCOM('Z012');
            
            std_file = 'CylGearPair 1 (spur gear).Z12';
%             file_name = sprintf('C:\\Program Files\\KISSsoft %s\\example\\%s', version, std_file);
            
            ex03.load_example_file(std_file);
            
            ex03.check_calculations();
            ex03.calculate();
            
            fprintf('Current profile shift gear 1 : %.3f\n', ex03.get_var('ZR[0].x.nul'));
            fprintf('Current tooth root safety gear 1: %.3f\n', ex03.get_var('ZPP[0].Fuss.SFnorm'));
            
            ex03.set_var('ZR[0].x.nul', 0.1);
            
            ex03.check_calculations();
            ex03.calculate();
            
            fprintf('New profile shift gear 1 : %.3f\n', ex03.get_var('ZR[0].x.nul'));
            fprintf('New tooth root safety gear 1: %.3f\n', ex03.get_var('ZPP[0].Fuss.SFnorm'));

        end
        
        function example_04()
            %EXAMPLE_04 from KUM
            %
            
            std_file = 'CylGearPair 1 (spur gear).Z12';
            ex04 = KISSsoftCOM('Z012');
            ex04.load_example_file(std_file);
            
            x = linspace(0.0, 1.0, 21);
            SF = zeros(size(x));
            
            for idx = 1:numel(x)
                ex04.set_var('ZR[0].x.nul', x(idx));
                SF(idx) = ex04.get_var('ZPP[0].Fuss.SFnorm');
            end
            
            figure;
            plot(x, SF, 'ko');
            xlabel('x, [-]');
            ylabel('S_F, [-]');
        end
        
    end
end