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
                obj.post   = actxserver(sprintf('Simpack.Post%s', default.version));
            else
                error(msg.identifier, "%s", msg.message);
            end
            
            if(strcmpi(default.mode, 'solver'))
                obj.COM = actxserver(sprintf('Simpack.SLV%s' , default.version));
            elseif(strcmpi(default.mode, 'GUI'))
                obj.COM = actxserver(sprintf('Simpack.GUI%s' , default.version));
            else
                error('SimpackCOM:mode', 'Valid modes are [SoLVer] and [GUI].');
            end
            
            file = dir(sprintf('@%s\\*.Nejad.spck', default.model_name));
            file_name = sprintf('%s\\%s', file.folder, file.name);
            obj.model = obj.open_model(file_name);
            
            obj.project = obj.post.Spck.addProject();
            
        end
        
        function mdl = open_model(obj, file)
            mdl = obj.COM.Spck.openModel(file);
        end
        
        function result = time_integration(obj, varargin)
            if(isempty(varargin))
                time_span = 1000.0;
            else
                time_span = varargin{1};
            end
            
            % Sets model time as 0
            obj.model.time.val = 0.0;
            obj.model.slv_active.val.integ_tend_duration.val = time_span;

%             tic;
            % time integration with measurements
            if(contains(obj.version, '2020'))
                result = obj.COM.Spck.Slv.integMeas(obj.current_model);
            elseif(contains(obj.version, '2018'))
                result = obj.COM.Spck.Slv.integ(obj.current_model);
                obj.COM.Spck.Slv.meas(obj.current_model, 1);
            end
%             disp(toc);
        end
        
        function delete(obj)
            obj.model.close();
            obj.post.Spck.closeProject(obj.project);
            obj.post.Spck.quit;
            obj.COM.Spck.Quit;
        end
        
%         function set_states_to_zero(obj)
%             Spck.currentModel.setStatesToZero(StateSet.JOINT_POS_STATE +
%             StateSet.JOINT_VEL_STATE +
%             StateSet.FLEXBODY_POS_STATE +
%             StateSet.FLEXBODY_VEL_STATE +
%             StateSet.CONSTR_FORCE_STATE +
%             StateSet.MARKER_ALG_STATE +
%             StateSet.FORCECONTR_ALG_STATE +
%             StateSet.FORCECONTR_DYN_STATE +
%             StateSet.FORCECONTR_DISC_STATE +
%             StateSet.FORCECONTR_DESCR_STATE +
%             StateSet.FORCECONTR_ROOT_STATE +
%             StateSet.BODY_POS_STATE +
%             StateSet.BODY_VEL_STATE +
%             StateSet.CONNECTION_POS_STATE +
%             StateSet.CONNECTION_VEL_STATE);
%         end
        
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
