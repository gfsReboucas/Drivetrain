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
            
            file = dir(sprintf('%s\\@NREL_5MW\\%s.spck', pwd, default.model_name));
            file_name = sprintf('%s\\%s', file.folder, file.name);
            obj.model = obj.open_model(file_name);
            
%             obj.project = obj.post.Spck.addProject();
            
        end
        
        function mdl = open_model(obj, file)
            mdl = obj.COM.Spck.openModel(file);
        end
        
        function save_file(obj)
            % delete problematic files first.
            path = fileparts(obj.model.origin);
            prob = dir(sprintf('%s/*.spck~*', path));
            
            for idx = 1:numel(prob)
                delete([prob(idx).folder, '\', prob(idx).name]);
            end
            
            obj.COM.Spck.currentModel.save();
        end
        
        function result = time_integration(obj)
              result = obj.COM.Spck.Slv.integMeas(obj.current_model);
%             % disabling translation and tilting for shafts:
%             svar = obj.find_element('$_switch_AIC');
%             svar.str.src = '1';
% 
%             obj.set_states_to_zero();
%             obj.apply_ICs();
%             
%             % enabling it again:
%             svar.str.src = '0';
%             obj.COM.Spck.Slv.equi(obj.current_model);
            
%             % time integration with measurements
%             if(contains(obj.version, '2018'))
%                 result = obj.COM.Spck.Slv.integ(obj.current_model);
%                 obj.COM.Spck.Slv.meas(obj.current_model, 1);
%             else
%                 result = obj.COM.Spck.Slv.integMeas(obj.current_model);
%             end
        end
        
        function [f_n, zeta] = modal_analysis(obj)
            obj.COM.Spck.Slv.eigen(obj.current_model, true);
            
            % Update state set:
            [path, name] = fileparts(obj.model.origin);
            load(sprintf('%s/output/%s.ev.mat', path, name), 'eigen');
            
            % resonances:
            f_n = double(eigen.eigenval.f0.values);
            zeta = double(eigen.eigenval.d.values);
%             f_nd = double(eigen.eigenval.f_d.values);

            [f_n, idx] = sort(f_n);
            zeta = zeta(idx);
            idx = f_n ~= 0.0;
            f_n = f_n(idx);
            zeta = zeta(idx);
        end
        
        function delete(obj)
            obj.model.close();
            obj.post.Spck.closeProject(obj.project);
            obj.post.Spck.quit;
            obj.COM.Spck.quit;
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
        
        function assemble(obj)
            obj.COM.Spck.Slv.assmbl(obj.current_model, true);
        end
        
        function elem = find_element(obj, name)
            elem = obj.current_model.findElement(name, true);
        end
        
        function val = get_subvar(obj, var_name)
            elem = obj.find_element(var_name);
            val = str2double(elem.str.src);
        end
        
        function set_subvar(obj, var_name, val)
            elem = obj.find_element(var_name);
            elem.str.src = num2str(val, 10);
        end
        
        function initial_step(obj)
            
%             fprintf('SCRIPT:\tSTART initial step on model [%s]...\n', upper(obj.current_model.name));
%             fprintf('SCRIPT:\t1. Disabling shaft translation and tilting...\n');
            obj.set_subvar('$SVG_loading.$_switch_AIC', 1);

            % Apply zero initial condtions:
%             fprintf('SCRIPT:\t2. Setting states to ZERO...\n');
            zero_IC = obj.find_element('$ST_all_zero');
            zero_IC.copyToModel();

            % Gear force elements: basic
%             fprintf('SCRIPT:\t2.1. Setting gear mesh force to type 204 (basic)...\n');
            obj.set_subvar('$SVG_loading.$_mesh_flag', 204);

%             fprintf('SCRIPT:\t2.2. Saving model...\n');
            obj.save_file();
            
            % Apply Initial Conditions solver:
%             fprintf('SCRIPT:\t3. Running Assemble System solver...\n');
            obj.assemble();

            % Update state set:
%             fprintf('SCRIPT:\t3.1. Saving model with updated state sets from solver...\n');
            state_AIC = obj.find_element('$ST_AIC');
            state_AIC.copyFromModel(false);
            obj.save_file();
            
            % Re-enabling translation and tilting for shafts:
%             fprintf('SCRIPT:\t4. Re-enabling shaft translation and tilting...\n');
            obj.set_subvar('$SVG_loading.$_switch_AIC', 0);

            % Calculating static equilibrium:
%             fprintf('SCRIPT:\t5. Calculating static equilibrium using:\n');
%             fprintf('SCRIPT:\t5.1. Time integration method...\n');
            solver = obj.find_element('$SLV_SolverSettings');
            solver.equi_st_meth.src = 1;
            solver.equi_cond.src = 0;
            obj.save_file();
            obj.COM.Spck.Slv.equi(obj.current_model, true, true);
            
            % Update state set:
            [path, name] = fileparts(obj.model.origin);
            full_name  = sprintf('%s/output/%s.spckst', path, name);
            
%             fprintf('SCRIPT:\t5.1.1. Updating state sets with result from solver (static equilibrium via time integration)...\n');
    
            static_eq = obj.find_element('$ST_static_eq');
            static_eq.importStateFile(full_name);
            static_eq.copyToModel();
            obj.save_file();
            
%             fprintf('SCRIPT:\t5.2. Newton method...\n');
            solver.equi_st_meth.src = 0;
            solver.equi_cond.src = 1;
            obj.save_file();
            obj.COM.Spck.Slv.equi(obj.current_model, true, true);

            % Update state set:
%             fprintf('SCRIPT:\t5.2.2. Updating state sets with results from solver (static equilibrium via Newton method and convergence w.r.t. acceleration)...\n');
            static_eq.importStateFile(full_name);
            static_eq.copyToModel();
            obj.save_file();
            
            % Gear force elements: basic
%             fprintf('SCRIPT:\t2.1. Setting gear mesh force to type 225 (advanced)...\n');
            obj.set_subvar('$SVG_loading.$_mesh_flag', 225);

%             fprintf('SCRIPT:\t6. Saving model...\n');
            obj.save_file();

%             fprintf('SCRIPT:\tEND initial_step.\n');
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
