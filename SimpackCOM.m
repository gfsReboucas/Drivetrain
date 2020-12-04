classdef SimpackCOM
    %SIMPACKCOM Wrapper class to deal with Simpack's COM interface. More
    % details about Simpack can be found on [1]. Run SimpackCOM.view_COM_methods() to
    % view all the methods provided by Simpack's COM interface. More info
    % about it at [2].
    %
    % [1] https://www.3ds.com/products-services/simulia/products/simpack/
    % [2] https://mathworks.com/help/matlab/ref/methodsview.html
    %
    %

% Apply state set to model:
% var obj1 = Spck.currentModel.findElement("$ST_all_zero");
% obj1.copyToModel();
% Spck.currentModel.save();

% Update state set from model:
% var obj2 = Spck.currentModel.findElement("$ST_AIC");
% obj2.copyFromModel();

    
    properties(Access = private)
        COM;
        post;
    end
    
    properties
        model;
        project;
    end
    
    properties(Dependent)
        current_model;
        version;
    end
    
    methods
        function obj = SimpackCOM(varargin)
            default = {'mode', 'SoLVer', ...
                       'version', '', ...
                       'model_name', 'NREL_5MW'};
                   
            [flag, msg] = SimpackCOM.is_installed();
            default = scaling_factor.process_varargin(default, varargin);
            
            if(flag)
                obj.post   = actxserver(sprintf('Simpack.Post.%s', default.version));
            else
                error(msg.identifier, "%s", msg.message);
            end
            
            if(strcmpi(default.mode, 'solver'))
                obj.COM = actxserver(sprintf('Simpack.Slv.%s' , default.version));
            elseif(strcmpi(default.mode, 'GUI'))
                obj.COM = actxserver(sprintf('Simpack.GUI.%s' , default.version));
            else
                error('SimpackCOM:mode', 'Valid modes are [SoLVer] and [GUI].');
            end
            
            file = dir(sprintf('@%s\\*.scaled.spck', default.model_name));
            file_name = sprintf('%s\\%s', file.folder, file.name);
            obj.model = obj.open_model(file_name);
            
            obj.project = obj.post.Spck.addProject();
            
        end
        
        function mdl = open_model(obj, file)
            mdl = obj.COM.Spck.openModel(file);
        end
        
        function result = time_integration(obj)
            
            % disabling translation and tilting for shafts:
            svar = obj.find_element('$_switch_AIC');
            svar.str.src = '1';

            obj.set_states_to_zero();
            obj.apply_ICs();
            
            % enabling it again:
            svar.str.src = '0';
            obj.COM.Spck.Slv.equi(obj.current_model);
            
            % time integration with measurements
            if(contains(obj.version, '2018'))
                result = obj.COM.Spck.Slv.integ(obj.current_model);
                obj.COM.Spck.Slv.meas(obj.current_model, 1);
            else
                result = obj.COM.Spck.Slv.integMeas(obj.current_model);
            end
        end
        
        function delete(obj)
            obj.model.close();
            obj.post.Spck.closeProject(obj.project);
            obj.post.Spck.quit;
            obj.COM.Spck.Quit;
        end
        
        function set_states_to_zero(obj)
            obj.model.setStatesToZero();
        end
        
        function apply_ICs(obj)
            obj.COM.Spck.Slv.applyInitialConditions(obj.current_model);
        end
        
        function pre_load(obj)
            obj.COM.Spck.Slv.preld(obj.current_model);
        end
        
        function elem = find_element(obj, name)
            elem = obj.current_model.findElement(name);
        end
        
    end
    
    methods(Static)
        function [val, msg] = is_installed()
            val = true;
            msg = [];
            try
                actxserver('Simpack.Post');
            catch msg
                val = false;
            end
        end
        
        function view_COM_methods()
            [flag, msg] = SimpackCOM.is_installed();
            if(flag)
                methodsview(actxserver('Simpack.Gui'));
                methodsview(actxserver('Simpack.Slv'));
                methodsview(actxserver('Simpack.Post'));
            else
                warning(msg.identifier, "%s", msg.message);
            end
        end
        
    end
    
    methods
        function val = get.current_model(obj)
            val = obj.COM.Spck.currentModel;
        end
        
        function val = get.version(obj)
            val = num2str(obj.COM.Spck.version.number);
        end
        
    end
    
end
