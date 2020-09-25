classdef (Abstract) Drivetrain
    %DRIVETRAIN This class implements SOME procedures for the dynamic
    % analysis and scaling of drivetrains. The safety factor for surface 
    % durability (pitting) is calculated according to ISO 6336 [1, 2]. The
    % NREL 5MW reference gearbox proposed by Nejad et. al. [4] is used as
    % the default drivetrain model, but other models should be implemented
    % in the future.
    %
    % References:
    % [1] ISO 6336-1:2006 Calculation of load capacity of spur and helical
    % gears -- Part 1: Basic principles, introduction and general influence
    % factors
    % [2] ISO 6336-2:2006 Calculation of load capacity of spur and helical
    % gears -- Part 2: Calculation of surface durability (pitting)
    % [3] IEC 61400-4:2012 Wind Turbines -- Part 4: Design Requirements for
    % wind turbine gearboxes.
    % [4] Nejad, A. R., Guo, Y., Gao, Z., Moan, T. (2016). Development of a
    % 5 MW reference gearbox for offshore wind turbines. Wind Energy. 
    % https://doi.org/10.1002/we.1884
    %
    % written by:
    % Geraldo Rebouças
    % - Geraldo.Reboucas@ntnu.no OR
    % - gfs.reboucas@gmail.com
    %
    % Postdoctoral Fellow at:
    %    Norwegian University of Science and Technology, NTNU
    %    Department of Marine Technology, IMT
    %    Marine System Dynamics and Vibration Lab, MD Lab
    %    https://www.ntnu.edu/imt/lab/md-lab
    %
    % see also NREL_5MW, DTU_10MW.
    
    properties(Access = public)
        stage         (1, :) Gear_Set;                                                                    % [class],  gearbox stages
        P_rated       (1, :)          {mustBeNumeric, mustBeFinite, mustBeNonnegative} = 5.0e3;              % [kW],     Rated power
        n_rotor       (1, :)          {mustBeNumeric, mustBeFinite, mustBeNonnegative} = 12.1;               % [1/min.], Rated rotor speed
        main_shaft    (1, :) Shaft                                                  = Shaft;              % [class],  Input Shaft
        m_Rotor       (1, :)          {mustBeNumeric, mustBeFinite, mustBeNonnegative} = 110.0e3;            % [kg],     Rotor mass according to [3]
        J_Rotor       (1, :)          {mustBeNumeric, mustBeFinite, mustBeNonnegative} = 57231535.0;         % [kg-m^2], Rotor mass moment of inertia according to [6]
        m_Gen         (1, :)          {mustBeNumeric, mustBeFinite, mustBeNonnegative} = 1900.0;             % [kg],     Generator mass according to [4]
        J_Gen         (1, :)          {mustBeNumeric, mustBeFinite, mustBeNonnegative} = 534.116;            % [kg-m^2], Generator mass moment of inertia [4]
        N_stage       (1, 1)          {mustBeInteger, mustBeFinite, mustBeNonnegative} = 3;                  % [-],      Number of stages
        dynamic_model (1, :)                                                        = @Dynamic_Formulation; % which dynamic model should be used to perform modal analysis on the Drivetrain.
        gamma                scaling_factor;                                                              % Scaling factors
        S_Hmin;      % [-], Minimum required safety factor for surface durability (pitting)
        S_Fmin;      % [-], Minimum required safety factor for tooth bending strength
        M; % Inertia matrix
        K; % Stiffness matrix
        D; % Damping matrix
        A; % State matrix
        load; % Load vector
        f_n; % natural frequencies
        f_nd; % damped natural frequencies
        mode_shape; % mode shapes
        mode_shape_d; % damped mode shapes
        zeta; % damping ratio
    end
    
    properties(Dependent)
        T_out;   % [N-m],    Output torque for each stage
        n_out;   % [1/min.], Output speed  for each stage
        u;       % [-],      Cumulative gear ratio
        S_H;     % [-],      Gear Safety factor for surface durability (against pitting)
        S_F;     % [-],      Gear Safety factor for bending strength
        S_shaft; % [-],      Shaft Safety factor against fatigue and yield
        S_G;     % [-],      Full safety factor for gear stages
    end
    
    properties(SetAccess = protected)
        % to store the values of some dependent variables:
%         S_H_val     (1, :) {mustBeNumeric, mustBeFinite, mustBePositive} = 1.25;   % [-],      Safety factor for surface durability (against pitting)
%         S_F_val     (1, :) {mustBeNumeric, mustBeFinite, mustBePositive} = 1.56;   % [-],      Safety factor for bending strength
        S_H_val     (1, :) {mustBeNumeric} = 1.25;   % [-],      Safety factor for surface durability (against pitting)
        S_F_val     (1, :) {mustBeNumeric} = 1.56;   % [-],      Safety factor for surface durability (against pitting)
        S_shaft_val (1, :) {mustBeNumeric, mustBeFinite, mustBePositive} = 1.0;    % [-],      Safey factor for the shafts
    end
    
    methods
        function obj = Drivetrain(varargin)
            
            % for default parameter:
            N_st = 3;
            stage = repmat(Gear_Set(), 1, N_st);
            for idx = 1:N_st
                stage(idx) = NREL_5MW.gear_stage(idx);
            end
            
            default = {'N_stage'      ,        3, ...
                       'P_rated'      ,        5.0e3, ...
                       'n_rotor'      ,       12.1, ...
                       'm_Rotor'      ,      110.0e3, ...
                       'J_Rotor'      , 57231535.0, ...
                       'm_Gen'        ,     1900.0, ...
                       'J_Gen'        ,      534.116, ...
                       'S_Hmin'       ,        1.25, ...
                       'S_Fmin'       ,        1.56, ...
                       'stage'        , stage, ...
                       'main_shaft'   , Shaft(), ...
                       'dynamic_model', @Dynamic_Formulation};
            
            default = scaling_factor.process_varargin(default, varargin);
            
            obj.P_rated    = default.P_rated;
            obj.n_rotor    = default.n_rotor;
            obj.main_shaft = default.main_shaft;
            
            obj.m_Rotor = default.m_Rotor;
            obj.J_Rotor = default.J_Rotor;
            
            obj.m_Gen = default.m_Gen;
            obj.J_Gen = default.J_Gen;
            
            if(isempty(default.stage))
                obj.N_stage = default.N_stage;
                for idx = 1:obj.N_stage
                    obj.stage(idx) = NREL_5MW.gear_stage(idx);
                end
            else
                obj.N_stage = length(default.stage);
                obj.stage = default.stage;
            end
            
            obj.S_Hmin = default.S_Hmin;
            obj.S_Fmin = default.S_Fmin;
            
            obj.dynamic_model = default.dynamic_model;
            dyn_calc = obj.dynamic_model(obj);
            obj.M = dyn_calc.M;
            obj.K = dyn_calc.K;
            obj.D = dyn_calc.D;
            obj.A = dyn_calc.A;
            obj.load = dyn_calc.c;
            
            [obj.f_n , obj.mode_shape  , ...
             obj.f_nd, obj.mode_shape_d, obj.zeta] = dyn_calc.modal_analysis();

        end
        
        function tab = disp(obj)
            %DISP display some properties of a Drivetrain object
            % description, symbol, unit, value
            if(isempty(obj))
                disp('\t0x0 empty Drivetrain object')
            else
                tab_str = {'Rated power',                                  'P',       'kW';
                           'Output Speed (Sun/Pinion)',                    'n_out',   '1/min.';
                           'Output Torque (Sun/Pinion)',                   'T_out',   'N-m';
                           'Minimum safety factor against pitting',        'S_Hmin',  '-';
                           'Safety factor against pitting (Sun/Pinion)',   'S_H1',    '-';
                           'Safety factor against pitting (Planet/Wheel)', 'S_H2',    '-';
                           'Safety factor (Shaft)',                        'S',       '-';
                           'Type',                                         '-',       '-';
                           'Normal module',                                'm_n',     'mm';
                           'Face width',                                   'b',       'mm';
                           'Center distance',                              'a_w',     'mm';
                           'Gear ratio',                                   'u',       '-';
                           'Number of planets',                            'p',       '-';
                           'Normal pressure angle',                        'alpha_n', 'deg.';
                           'Helix angle',                                  'beta',    'deg.';
                           'Number of teeth (Sun/Pinion)',                 'z_1',     '-';
                           'Number of teeth (Planet/Wheel)',               'z_2',     '-';
                           'Number of teeth (Ring)',                       'z_3',     '-';
                           'Profile shift coefficient (Sun/Pinion)',       'x_1',     '-';
                           'Profile shift coefficient (Planet/Wheel)',     'x_2',     '-';
                           'Profile shift coefficient (Ring)',             'x_3',     '-';
                           'Reference diameter (Sun/Pinion)',              'd_1',     'mm';
                           'Reference diameter (Planet/Wheel)',            'd_2',     'mm';
                           'Reference diameter (Ring)',                    'd_3',     'mm';
                           'Mass (Sun/Pinion)',                            'm_1',     'kg';
                           'Mass (Planet/Wheel)',                          'm_2',     'kg';
                           'Mass (Ring)',                                  'm_3',     'kg';
                           'Mass mom. inertia (Sun/Pinion)',               'J_xx1',   'kg-m^2';
                           'Mass mom. inertia (Planet/Wheel)',             'J_xx2',   'kg-m^2';
                           'Mass mom. inertia (Ring)',                     'J_xx3',   'kg-m^2';
                           'Mass mom. inertia (Sun/Pinion)',               'J_yy1',   'kg-m^2';
                           'Mass mom. inertia (Planet/Wheel)',             'J_yy2',   'kg-m^2';
                           'Mass mom. inertia (Ring)',                     'J_yy3',   'kg-m^2';
                           'Mass mom. inertia (Sun/Pinion)',               'J_zz1',   'kg-m^2';
                           'Mass mom. inertia (Planet/Wheel)',             'J_zz2',   'kg-m^2';
                           'Mass mom. inertia (Ring)',                     'J_zz3',   'kg-m^2';
                           };
                
                Parameter = tab_str(:, 1);
                Symbol    = tab_str(:, 2);
                Unit      = tab_str(:, 3);
                
                tab_val = cell(obj.N_stage + 2, 1);
                
                tab_val{1} = table(Parameter, Symbol);
                
                for idx = 1:obj.N_stage
                    val_stg = {obj.P_rated;
                               obj.n_out(idx);
                               obj.T_out(idx);
                               obj.S_Hmin;
                               obj.S_H(2*idx - 1);
                               obj.S_H(2*idx);
                               obj.S_shaft(idx + 1);
                               obj.stage(idx).configuration;
                               obj.stage(idx).m_n;
                               obj.stage(idx).b;
                               obj.stage(idx).a_w;
                               obj.stage(idx).u;
                               obj.stage(idx).N_p;
                               obj.stage(idx).alpha_n;
                               obj.stage(idx).beta;
                               obj.stage(idx).z(1);
                               obj.stage(idx).z(2);
                               '-*-';
                               obj.stage(idx).x(1);
                               obj.stage(idx).x(2);
                               '-*-';
                               obj.stage(idx).d(1);
                               obj.stage(idx).d(2);
                               '-*-';
                               obj.stage(idx).mass(1);
                               obj.stage(idx).mass(2);
                               '-*-';
                               obj.stage(idx).J_x(1);
                               obj.stage(idx).J_x(2);
                               '-*-';
                               obj.stage(idx).J_y(1);
                               obj.stage(idx).J_y(2);
                               '-*-';
                               obj.stage(idx).J_z(1);
                               obj.stage(idx).J_z(2);
                               '-*-';
                               };
                    
                    if(strcmp(obj.stage(idx).configuration, 'planetary'))
                        val_stg{18} = obj.stage(idx).z(3);
                        val_stg{21} = obj.stage(idx).x(3);
                        val_stg{24} = obj.stage(idx).d(3);
                        val_stg{27} = obj.stage(idx).mass(3);
                        val_stg{30} = obj.stage(idx).J_x(3);
                        val_stg{33} = obj.stage(idx).J_y(3);
                        val_stg{36} = obj.stage(idx).J_z(3);
                    end
                    
                    tab_val{idx + 1} =  table(val_stg, 'variableNames', sprintf("Stage_%d", idx));
                end
                
                tab_val{idx + 2} = table(Unit);
                
                tab = [tab_val{:}];
                
                if(nargout == 0)
                    fprintf('Gear stages:\n');
                    disp(tab);
                    fprintf('Main shaft:\n');
                    disp(obj.main_shaft);
                    clear tab;
                end
            end            
        end
        
        function tab = comparison(ref, sca)
            tab_stg = cell(ref.N_stage, 1);
            
            for idx = 1:(ref.N_stage)
                [~, tab_str] = stage_comparison(ref, sca, idx);
                tab_tmp = table(tab_str(:, 4), tab_str(:, 5), tab_str(:, 6), 'variableNames', ['Reference', 'Scale', 'Ratio']);
                tab_stg{idx} = table(tab_tmp, 'variableNames', sprintf('Stage_%d', idx));
            end
            
            tab_left =  table(tab_str(:, 1), tab_str(:, 2), 'variableNames', ['Parameter', 'Symbol']);
            tab_right = table(tab_str(:, 3), 'variableNames', 'Unit');
            
            tab = [tab_left, tab_stg{:}, tab_right];
            
            if(nargout == 0)
                disp(tab);
                clear tab;
            end
            
        end
        
        function [tab, tab_str] = stage_comparison(ref, sca, idx)
            if(ref.N_stage ~= sca.N_stage)
                error('Drivetrain:different_stages', 'different number of gear stages.');
            elseif(idx > ref.N_stage)
                error('Drivetrain:big_idx', 'idx is bigger than the number of stages.');
            elseif(idx < 1)
                error('Drivetrain:negative_idx', 'idx is smaller than 1.');
            end
            
            tab_str = {'Rated power',                                  'P',       'kW',     ref.P_rated,                sca.P_rated,                sca.P_rated               /ref.P_rated;                % 1
                       'Output Speed (Sun/Pinion)',                    'n_out',   '1/min.', ref.n_out(idx),             sca.n_out(idx),             sca.n_out(idx)            /ref.n_out(idx);             % 2
                       'Output Torque (Sun/Pinion)',                   'T_out',   'N-m',    ref.T_out(idx),             sca.T_out(idx),             sca.T_out(idx)            /ref.T_out(idx);             % 3
                       'Safety factor against pitting (Sun/Pinion)',   'S_H1',    '-',      ref.S_H(2*idx - 1),         sca.S_H(2*idx - 1),         sca.S_H(2*idx - 1)        /ref.S_H(2*idx - 1);         % 4
                       'Safety factor against pitting (Planet/Wheel)', 'S_H2',    '-',      ref.S_H(2*idx),             sca.S_H(2*idx),             sca.S_H(2*idx)            /ref.S_H(2*idx);             % 5
                       'Safety factor (Shaft)',                        'S',       '-',      ref.S_shaft(idx + 1),       sca.S_shaft(idx + 1),       sca.S_shaft(idx + 1)      /ref.S_shaft(idx + 1);
                       };
            
            [~, str_stg] = comparison(ref.stage(idx), sca.stage(idx));
            tab_str = [tab_str; str_stg];
            
            Parameter = tab_str(:, 1);
            Symbol    = tab_str(:, 2);
            Unit      = tab_str(:, 3);
            Reference = tab_str(:, 4);
            Scale     = tab_str(:, 5);
            Ratio     = tab_str(:, 6);
            
            tab = table(Parameter, Symbol, Scale, Reference, Ratio, Unit);
            
            if(nargout == 0)
                disp(tab);
                clear tab tab_str;
            end
            
        end
        
        function plot(obj)
            hold on;
            axis equal;
            box on;
            for idx = 1:obj.N_stage
                subplot(1, obj.N_stage, idx)
                obj.stage(idx).plot;
            end
        end
        
        function rectangle(obj, varargin)
            if(nargin == 1)
                C_0 = zeros(2, 1);
            else
                C_0 = varargin{1};
            end
            
            % LINSPECER: Plot lots of lines with very distinguishable and 
            % aesthetically pleasing colors. It can be dowloaded from
            % MATLAB's File Exchange on:
            % https://se.mathworks.com/matlabcentral/fileexchange/42673-beautiful-and-distinguishable-line-colors-+-colormap
            color = linspecer(6, 'qualitative');
            
            hold on;
            C_s = [obj.main_shaft.L/2.0 0.0]' + C_0;
            rectangle(obj.main_shaft, C_s, color(6, :), 'edgeColor', 'k', 'lineStyle', '-' , 'faceColor', color(6, :));
            
            C_x = [obj.main_shaft.L + C_0(1);
                   obj.stage(1).carrier.b + obj.stage(1).output_shaft.L;
                   obj.stage(2).carrier.b + obj.stage(2).output_shaft.L];

            for idx = 1:obj.N_stage
                rectangle(obj.stage(idx), [sum(C_x(1:idx)) C_0(2)]');
            end
            hold off;
        end
        
        function plot_comp(DT1, DT2)
            if(DT1.N_stage ~= DT2.N_stage)
                error('Drivetrain:different_stage_number', 'Both DT''s should have the same number of stages.');
            end
            
            hold on;
            axis equal;
            box on;
            for idx = 1:DT1.N_stage
                subplot(2, DT1.N_stage, idx)
                DT1.stage(idx).plot;
                title(sprintf('Stage %d', idx));
                subplot(2, DT2.N_stage, idx + DT2.N_stage)
                DT2.stage(idx).plot;
            end

        end
        
    end
    
    %% Calculations:
    methods
        %% Dynamics:
        function [fn, modeShape] = resonances(obj, varargin)
            %RESONANCES returns the first N resonances and mode shapes of a
            % Drivetrain object. The resonances can be normalized or not.
            %

            default = {'N'        , calc.n_DOF(end), ...
                       'normalize', false};

            default = scaling_factor.process_varargin(default, varargin);
            
            N = default.N;
            normalize = default.normalize;
            
            fn = obj.f_n;
            modeShape = obj.mode_shape;
            
            N_fn = numel(fn);
            
            if(N <= 0)
                error('N = %d < 0. It should be positive and smaller than %d.', N, N_fn);
            elseif(N > N_fn)
                error('N = %d > %d. It should be positive and smaller than %d.', N, N_fn, N_fn);
            else
                fn = fn(1:N);
                modeShape = modeShape(:, 1:N);
            end
            
            if(normalize == true)
                fn = fn./fn(1);
                fn = fn(2:end);
            end
            
        end
        
        function [f, modeShape] = nth_resonance(obj, n)
            [fn, modeShape] = obj.resonances('N', n, 'normalize', false);
            
            if((n < 1) || (n > numel(fn)))
                error('n = %d > or ~= 0.');
            end
            
            f = fn(n);
            modeShape = modeShape(:, n);
        end
        
        function [model, result] = Simpack_time_integration(obj)
            sim = SimpackCOM();
            
%             file = dir(sprintf('@%s\\*.Nejad.spck', class(obj)));
%             file_name = sprintf('%s\\%s', file.folder, file.name);
%             model = sim.open_model(file_name);
            model = sim.model;
            result = sim.time_integration();
        end
        
        function SIMPACK_version(obj, varargin)
            if(nargin == 1)
                version = '2018';
                task = 'lsa';
            elseif(nargin == 2)
                version = '2018';
                task = varargin{1};
            elseif(nargin == 3)
                version = varargin{1};
                task = varargin{2};
            else
                error('too many arguments.');
            end
            
            % create SIMPACK COM solver interface:
            solver = actxserver(sprintf('Simpack.Slv.%s', version));
            
            % Open model:
            file = dir(sprintf('@%s\\*.spck', class(obj)));
            file_name = sprintf('%s\\%s', file.folder, file.name);
            model = solver.Spck.openModel(file_name);
            
            switch(task)
                case 'lsa'
                    % Linear System Analysis:
                    solver.Spck.Slv.lsa(model);
                    
                    % the following part seems useless:
                    % create SIMPACK COM result interface:
                    post = actxserver(sprintf('Simpack.Post.%s', version));

                    % add project:
                    project = post.Spck.addProject();
                    file = dir(sprintf('@%s\\*\\*.lsa.sbr', class(obj)));
                    file_name = sprintf('%s\\%s', file.folder, file.name);

                    % add Linear System Response result file to project:
                    result = project.addResultFile(file_name);

                case 'eigen'
                    % Modal Analysis:
                    result = solver.Spck.Slv.eigen(model, false);
                    nx = result.numEigenvalues; % number of eigenvalues
                    freq = zeros(nx, 1);

                    for idx = 1:nx
                        freq(idx) = result.freq(idx - 1);
                    end

                    idx = find(freq > 1.0e-5, 1, 'first');
                    freq = freq(idx:end);
                    
                case 'ssm'
                    % State Space Matrix:
                    result = solver.Spck.Slv.ssm(model, ...
                                                 2, ... % 'new' MATLAB format
                                                 false); % don't re-use an existing solver
                    nx = result.stateDim;
                    nu = result.inputDim;
                    ny = result.outputDim;
                    
                    A = zeros(nx, nx);        B = zeros(nx, nu);
                    C = zeros(ny, nx);        D = zeros(ny, nu);
                    
                    n = max([nx nu ny]);
                    
                    for row = 1:n
                        for col = 1:n
                            if((row <= nx) && (col <= nx))
                                A(row, col) = result.A(row - 1, col - 1);
                            elseif((row <= nx) && (col <= nu))
                                B(row, col) = result.B(row - 1, col - 1);
                            elseif((row <= ny) && (col <= nx))
                                C(row, col) = result.C(row - 1, col - 1);
                            elseif((row <= ny) && (col <= nu))
                                D(row, col) = result.D(row - 1, col - 1);
                            end
                        end
                    end
            end
            
        end
        
        %% Pitting:
        [SH, SF, SShaft, calc] = safety_factor_stage(obj, idx)
        
        function [SH_vec, SF_vec, SShaft_vec, calc] = safety_factors(obj)
            %SAFETY_FACTORS returns the safety factors of the shafts and
            % gear stages from a Drivetrain object. The shafts' safety
            % factor is calculated against fatigue and yield according to
            % Shigley's book and the gears safety factor is calculated
            % against pitting and teeth root failure according to 
            % ISO 6336.
            %
            
            K_f  = 1.0; % [-], Fatigue stress-concentration factor for bending
            K_fs = 1.0; % [-], Fatigue stress-concentration factor for torsion
            
            T_m  = obj.T_out(1)*obj.stage(1).u;

            max_Np = max([obj.stage.N_p]) + 1;
            
            SH_vec     = zeros(obj.N_stage    , max_Np);
            SF_vec     = zeros(obj.N_stage    , max_Np);
            SShaft_vec = zeros(obj.N_stage + 1, 1);
            
            SShaft_vec(1) = obj.main_shaft.safety_factors(K_f, K_fs, T_m);

            jdx = 0;
            for idx = 1:obj.N_stage
                [SH, SF, SShaft_vec(idx + 1), calc] = obj.safety_factor_stage(idx);
                
%                 calc.save_KS(sprintf('%s\\@%s\\test_0%d', pwd, class(obj), idx));
                clear('calc');

                kdx = jdx + (1:length(SH));
                jdx = kdx(end);
                
                SH_vec(kdx) = SH;
                SF_vec(kdx) = SF;
            end
            
            SH_vec(SH_vec == 0) = [];
            SF_vec(SF_vec == 0) = [];
        end
        
        %% Scaling:
        function [gamma_opt, res, stage_sca] = scaled_version_stage(obj, idx, gamma_idx, aspect, gamma_P, gamma_n)
            [SHref, SFref] = obj.safety_factor_stage(idx);
            SGref = [SHref, SFref];
            SGref(SGref == 0) = [];
            SGref(isnan(SGref)) = [];
            
            Smin = [~isnan(SHref)*obj.S_Hmin, ~isnan(SFref)*obj.S_Fmin];
            Smin = Smin(Smin ~= 0.0);
            obj_tmp = obj;
            obj_tmp.P_rated = gamma_P*obj.P_rated;
            obj_tmp.n_rotor = gamma_n*obj.n_rotor;
            
            vec = [gamma_idx('m_n') gamma_idx('b')]';
            if(strcmp(aspect, 'integrity'))
                gamma_0 = mean(vec);
            elseif(strcmp(aspect, 'integrity_refine'))
                gamma_0 = vec;
            end
            
            gamma_Max = ones(size(gamma_0));
            gamma_min = gamma_Max*1.0e-6;
            
            fun_asp = @(x)obj_tmp.scaled_safety_factor_stage(x, idx, SGref);
            fun_min = @(x)(norm(fun_asp(x))^2);
            
            fun_ineq = @(x)(Smin - obj_tmp.scaled_safety_factor_stage(x, idx, SGref*0));
            
            constraint_fun = @(x)deal(fun_ineq(x), fun_asp(x)); % inequalities, equalities
            
            opt_solver = optimoptions('fmincon', 'display', 'notify');
            
            [gamma_opt, res] = fmincon(fun_min, gamma_0, [], [], [], [], ...
                gamma_min, gamma_Max, constraint_fun, opt_solver);
            
            if(strcmp(aspect, 'integrity'))
                gm = scaling_factor({'m_n', gamma_opt, ...
                                     'b'  , gamma_opt});
            elseif(strcmp(aspect, 'integrity_refine'))
                gm = scaling_factor({'m_n', gamma_opt(1), ...
                                     'b'  , gamma_opt(2)});
            end
            
            stage_sca = obj.stage(idx).scale_by(gm);
        end
        
        function S_diff = scaled_safety_factor_stage(obj, gamma_val, idx, SGref)
            gm = scaling_factor({'m_n', 1.0, ...
                                 'b'  , 1.0});
            
            if(length(gamma_val) == 1)
                gm('m_n') = gamma_val;
                gm('b')   = gamma_val;
            elseif(length(gamma_val) == 2)
                gm('m_n') = gamma_val(1);
                gm('b')   = gamma_val(2);
            end
                
            obj.stage(idx) = obj.stage(idx).scale_by(gm);
            
            [SH, SF] = obj.safety_factor_stage(idx);
            SG = [SH SF];
            SG(SG == 0) = [];
            SG(isnan(SG)) = [];
            
            S_diff = SG - SGref;
        end
        
        function obj_sca = scale_all(obj_ref, gamma_P, gamma_n, gamma)
            %SCALE_ALL returns an object scaled by the factors gamma_P and 
            % gamma_n for its rated power and rotor speed and gamma for 
            % normal module, face width, shaft dimensions (diameter and 
            % length), mass and mass moment of inertia of rotor and 
            % generator.
            %
            % current limitations/constraints:
            % - shaft diameter is proportional to the cubic root of torque;
            %
            
            if((length(gamma) ~= length(obj_ref.gamma)) || (~isa(gamma, 'scaling_factor')))
                error('gamma must be a scaling_factor object with %d elements.', length(obj_ref.gamma));
            end
            
            key_set = obj_ref.gamma.name;
            
            gamma_mR  = 1.0; %gamma('m_R');
            gamma_JR  = gamma('J_R');
            gamma_mG  = 1.0; %gamma('m_G');
            gamma_JG  = gamma('J_G');
%             gamma_Ls  = gamma('L_S');
%             gamma_ds  = gamma('d_S');
            
            gamma_Torque = gamma_P/gamma_n; % Applied torque
            
            % Main shaft:
%             gamma_Ls = gamma('L_3'); % length
            gamma_ds = nthroot(gamma_Torque, 3.0);
            
            % Shaft diameters are scaled in the same way
            idx_dd = contains(key_set, 'd');
            idx_LL = contains(key_set, 'L');
            
            gamma(key_set(idx_dd)) = gamma_ds*ones(sum(idx_dd), 1);
            
            gamma_dd = gamma(key_set(idx_dd));
            gamma_LL = gamma(key_set(idx_LL));
            
            stage_sca = [Gear_Set, Gear_Set, Gear_Set];
            
            for idx = 1:obj_ref.N_stage
                jdx = 4*idx + (-3:0);
                gamma_stage = gamma(key_set(jdx));
                
                stage_sca(idx) = obj_ref.stage(idx).scale_aspect(gamma_stage, 'Gear_Set');
            end
            
            idx = idx + 1;
            
            main_shaft_sca = Shaft(obj_ref.main_shaft.d*gamma_dd(idx), ...
                                   obj_ref.main_shaft.L*gamma_LL(idx));
            
            obj_sca = Drivetrain(obj_ref.N_stage, ...
                                 stage_sca, ...
                                 obj_ref.P_rated *gamma_P, ...
                                 obj_ref.n_rotor *gamma_n, ...
                                 main_shaft_sca, ...
                                 obj_ref.m_Rotor*gamma_mR, ...
                                 obj_ref.J_Rotor*gamma_JR, ...
                                 obj_ref.m_Gen  *gamma_mG, ...
                                 obj_ref.J_Gen  *gamma_JG);
                             
            obj_sca.dynamic_model = obj_ref.dynamic_model;
        end
        
        function obj_sca = scale_aspect(obj_ref, gm_val, gamma_T, aspect)
            %SCALE_ASPECT returns an object scaled by the factors gamma_P 
            % and gamma_n for its rated power and rotor speed. The
            % argument gamma can be used to scale the Drivetrain object
            % considering various aspects.
            %
            
            key_set = obj_ref.gamma.name;
            
            gamma_full = obj_ref.gamma;
            
            asp_tmp = reshape(aspect, 1, numel(aspect));
            flag = any(ismember(key_set, asp_tmp));
            
            if(~flag)
                error('Aspect contains NO parameter related to %s.', class(obj_ref));
            end
            
            lin = size(aspect, 1);
            
            if(numel(gm_val) ~= lin)
                error('gm_val should have the same number of lines of aspect.')
            end
            
            for idx = 1:lin
                asp_tmp = aspect(idx, :);
                asp_tmp = asp_tmp(asp_tmp ~= '*');
                col = numel(asp_tmp);
                
                for jdx = 1:col
                    gamma_full(aspect(idx, jdx)) = gm_val(idx, :);
                end
            end
            
            % Scaling the shaft's diameter:
            gamma_ds = nthroot(gamma_T, 3.0);

            idx_D = contains(key_set, 'd');
            gamma_full(key_set(idx_D)) = gamma_ds*ones(sum(idx_D), 1);
            
            if(isa(obj_ref, 'Drivetrain'))
                obj_sca = scale_all(gamma_P, gamma_n, gamma_full);
            else
                scale_func = str2func(class(obj_ref));
                obj_sca = scale_func(gamma_P, gamma_n, gamma_full);
            end
            
            % Examples of aspects:
            % stage:
            % ["m_n1", "b_1"; ...
            %  "m_n2", "b_2"; ...
            %  "m_n3", "b_3"];
            % gear:
            % reshape(stage', numel(stage), 1);
            % K_MMI:
            % ["L_1", "L_2", "L_3", "L_S"; ...
            %  "J_R", "J_G", "*",   "*" ]; % "*" are used to complete the
            %  string array.
            %
            
            % [TODO]:
            % 'K_MMI_stage': scales parameters related to the stiffness and
            % mass moment of inertia, including the gear stages:
            % min |1 - f(sca)/f(ref)|^2 = f(x)
            % subjected to: S(min) - S(sca) <= 0 = g(x)
            %               S(ref) - S(sca)  = 0 = h(x)
            % ["m_n1", "b_1"; ...
            %  "m_n2", "b_2"; ...
            %  "m_n3", "b_3"; ...
            %  "L_1" , "L_2", "L_3"; ...
            %  "J_R", "J_G"; ...
            %  "L_S"];
            % 'K_MMI_gear'
            % 'M_MMI_stage_det'
            % 'K_MMI_gear_det'
            % 'K_shaft' : scales only parameters related to the shaft's 
            % 'K': scales only stiffness-related parameters. The mesh
            % stiffnesss is assumed to be proportinal to the gear's face
            % width.

        end
        
        function [obj_sca, gamma_opt, res, gamma_sep] = scaled_version(obj_ref, P_scale, varargin)
            %SCALED_VERSION returns an scaled object together with the 
            % scaling factor for its parameters and residual error from the
            % scaling optimization process. Different scaling can be
            % obtained by using:
            % - aspect_set: Different model approaches;
            % - normalize_freq: Normalizing or not the target resonances;
            % - N_freq: choosing how many resonances should be considered;
            % - gamma_0: initial value of gamma;
            % during the scaling optimization process;
            % Additionaly, one can define an initial value for gamma to be
            % used on the optimization process
            
            gamma_P = P_scale/obj_ref.P_rated;
            gm0 = scaling_factor(obj_ref.gamma.name, ones(length(obj_ref.gamma), 1)*nthroot(gamma_P, 3));
            
            constraint = ["L_1", "L_2";
                          "L_3", "L_S";
                          "J_R", "J_G"];
            
            default = {'n_scale'       , obj_ref.n_rotor, ...
                       'normalize_freq', true, ...
                       'N_freq'        , 10, ...
                       'aspect_set'    , 'DA', ...
                       'gamma_0'       , gm0, ...
                       'constraint'    , constraint};
            
            default = scaling_factor.process_varargin(default, varargin);
            
            constraint = default.constraint;
            
            gamma_0 = default.gamma_0;
            
            % Input scaling factors:
            gamma_0('P') = gamma_P;
            gamma_0('n') = default.n_scale/obj_ref.n_rotor;
            gamma_T = gamma_0('P')/gamma_0('n');
            
            % Scaling factors from dimensional analysis:
            gamma_length = nthroot(gamma_T,     3.0);
            
            gamma_MMI    = power(gamma_0('P'), 5.0/3.0);
            
            fprintf('Scaling the %s drivetrain to %.2f %% of its reference rated power (= %.1f kW).\n', class(obj_ref), gamma_0('P')*[100.0, obj_ref.P_rated]);
            
            id = {'MATISO_6336:KV', ...
                  'MATISO_6336:SF', ...
                  'MATISO_6336:KS', ...
                  'MATLAB:nearlySingularMatrix', ...
                  'MATLAB:COM:InvalidProgid', ...
                  'Dynamic_Formulation:imag', ...
                  'Dynamic_Formulation:RB', ...
                  'CVX:Renamed'};
            
            for idx = 1:length(id)
                warning('off', id{idx});
            end
            opt_wrt = 'Optimizing w. r. t. [%s]...\n';
            
            %% 1. Stage scaling:
            aspect_set = default.aspect_set;
            aspect_1 = 'integrity';
            
            if(any(aspect_set == aspect_1))
                fprintf(opt_wrt, upper(aspect_1));
                
                gm_val_1 = zeros(obj_ref.N_stage, 1);
                res_1    = zeros(obj_ref.N_stage, 1);
                
                for idx = 1:obj_ref.N_stage
                    gamma_idx = gamma_0.ends_with(num2str(idx));

                    [gm_val_1(idx, :), res_1(idx, :)] = ...
                        obj_ref.scaled_version_stage(idx, gamma_idx, aspect_1, gamma_0('P'), ...
                                                                               gamma_0('n'));
                end
                
            elseif(any(aspect_set == 'DA'))
                m_n_tmp =[obj_ref.stage.m_n]';
                gm_val_1 = Rack.module(m_n_tmp*gamma_length, 'calc', 'nearest')./m_n_tmp;
                
                res_1 = Inf(3, 1);
            else
                gm_val_1 = [mean(gamma_0(["m_n1" "b_1"]));
                            mean(gamma_0(["m_n2" "b_2"]));
                            mean(gamma_0(["m_n3" "b_3"]))];
                       
                res_1 = Inf(3, 1);
            end
            
            for idx = 1:obj_ref.N_stage
                gamma_0(sprintf('m_n%d', idx)) = gm_val_1(idx);
                gamma_0(sprintf('b_%d' , idx)) = gm_val_1(idx);
            end
            
            gm_val_1  = reshape(repmat(gm_val_1, 1, 2)', 2*length(gm_val_1), 1);
            res_1     = reshape(repmat(res_1   , 1, 2)', 2*length(res_1)   , 1);
            
            %% 2. Gears:
            aspect_2 = 'integrity_refine';
            
            if(any(aspect_set == aspect_2))
                fprintf(opt_wrt, upper(aspect_2));
                
                gm_val_2 = zeros(2*obj_ref.N_stage, 1);
                res_2    = zeros(2*obj_ref.N_stage, 1);
                
                for idx = 1:obj_ref.N_stage
                    gamma_idx = gamma_0.ends_with(num2str(idx));
                    
                    jdx = 2*idx + (-1:0);
                    [gm_val_2(jdx, :), res_2(jdx, :)] = ...
                        obj_ref.scaled_version_stage(idx, gamma_idx, aspect_2, gamma_0('P'), ...
                                                                               gamma_0('n'));
                end

            elseif(any(aspect_set == 'DA'))
                m_n_tmp = [obj_ref.stage.m_n]';
                gamma_mn = Rack.module(m_n_tmp*gamma_length, 'calc', 'nearest')./m_n_tmp;
                
                gm_val_2          = ones(6, 1)*gamma_length;
                gm_val_2(1:2:end) = gamma_mn;
                
                res_2 = Inf(6, 1);
            else
                gm_val_2 = gamma_0(["m_n1" "b_1" "m_n2" "b_2" "m_n3" "b_3"]);
                res_2 = Inf(size(gm_val_2));
            end
            
            if(all(isinf([res_1' res_2'])))
                gm_val_12 = gm_val_1;
                res_12 = res_1;
            else
                idx_min = res_1 <= res_2;

                gm_val_12 = diag(idx_min)*gm_val_1 + diag(~idx_min)*gm_val_2;
                res_12   = diag(idx_min)*res_1   + diag(~idx_min)*res_2;
            end
            
            for idx = 1:obj_ref.N_stage
                jdx = 2*idx + (-1:0);
                q = jdx(1);
                w = jdx(end);
                gamma_0(sprintf('m_n%d', idx)) = gm_val_12(q);
                gamma_0(sprintf('b_%d' , idx)) = gm_val_12(w);
            end
            
            par_d = gamma_0.name_contains('d_');
            
            gamma_0(par_d) = gamma_length;
            gm0_cell = gamma_0.to_cell();
            
            scale_func = str2func(class(obj_ref));
            obj_12 = scale_func(gm0_cell{:});
            
            %% 3. Shaft stiffness and Mass moment of inertia of rotor and generator:
            f_n_ref = obj_ref.resonances('N', default.N_freq, 'normalize', default.normalize_freq);
            
            aspect_3 = 'resonances';
            
            opt_solver = optimoptions('fmincon', 'display', 'notify');

            par_L = obj_12.gamma.name_contains('L_')';
            par_J = obj_12.gamma.name_contains('J_')';
            
            if(any(aspect_set == aspect_3))
                fprintf(opt_wrt, upper(aspect_3));
                
                param_constraint = [par_L;
                                    par_J, repmat("*", 1, abs(length(par_L) - length(par_J)))];
                                
                fun_asp  = @(x)(1.0 - obj_12.scale_resonances(x, gamma_0, ...
                                                                 default.N_freq, ...
                                                                 default.normalize_freq, ...
                                                                 param_constraint)./f_n_ref);
                fun_min = @(x)(norm(fun_asp(x))^2);
                
                constraint_fun = @(x)deal([], fun_asp(x)); % inequalities, equalities

                gm_03 = [mean(gamma_0(par_L)); ... % length
                         mean(gamma_0(par_J))]; % mass mom. inertia

                gamma_Max = ones(size(gm_03));
                gamma_min = gamma_Max*1.0e-6;

                [gm_val_3, res_3, ~] = fmincon(fun_min, gm_03, [], [], [], [], ...
                                               gamma_min, gamma_Max, ...
                                               constraint_fun, opt_solver);
            elseif(any(aspect_set == 'DA'))
                gm_val_3 = [gamma_length;
                            gamma_MMI];
                res_3 = Inf;
            else
                gm_val_3 = [mean(gamma_0(["L_1" "L_2" "L_3" "L_S"]));
                            mean(gamma_0(["J_R" "J_G"]))];
                       
                res_3 = Inf;
            end
            
            gamma_0(par_L) = gm_val_3(1)*ones(length(par_L), 1);
            gamma_0(par_J) = gm_val_3(2)*ones(length(par_J), 1);

            %% 4. Detailed version of 3:
            aspect_4 = 'resonances_refine';
            
            if(any(aspect_set == aspect_4))
                fprintf(opt_wrt, upper(aspect_4));
                
                % function for aspect_4: gamma_P and gamma_n are set to 1.0
                % because the drivetrain is already scaled for these
                % parameters.
                fun_asp  = @(x)(1.0 - obj_12.scale_resonances(x, gamma_0, ...
                                                                 default.N_freq, ...
                                                                 default.normalize_freq, ...
                                                                 constraint)./f_n_ref);
                fun_min = @(x)(norm(fun_asp(x))^2);

                constraint_fun = @(x)deal([], fun_asp(x)); % inequalities, equalities

                n = size(constraint, 1);
                
                gm_04 = zeros(n, 1);
                
                for idx = 1:n
                    gm_04(idx) = mean(gamma_0(constraint(idx,:)));
                end
                
                gamma_Max = ones(n, 1);
                gamma_min = gamma_Max*1.0e-6;

                [gm_val_4, res_4, ~] = fmincon(fun_min, gm_04, [], [], [], [], ...
                                               gamma_min, gamma_Max, ...
                                               constraint_fun, opt_solver);
            elseif(any(aspect_set == 'DA'))
                gm_val_4      = gamma_length*ones(6, 1);
                gm_val_4(4:5) = gamma_MMI; 
                
                res_4 = Inf;
            else
                gm_val_4 = gamma_0(["L_1" "L_2" "L_3" "L_S" "J_R" "J_G"]);
                res_4    = Inf;
            end
            
            if(all(isinf([res_3 res_4])) || (res_3 <= res_4))
                gm_val_34 = gm_val_3;
                res_34 = res_3;
                
                gamma_0(par_L) = gm_val_34(1)*ones(length(par_L), 1);
                gamma_0(par_J) = gm_val_34(2)*ones(length(par_J), 1);
            else
                gm_val_34 = gm_val_4;
                res_34 = res_4;
                
                for idx = 1:n
                    asp_tmp = constraint(idx, :);
                    asp_tmp(asp_tmp == '*') = [];
                    
                    m = size(asp_tmp);
                    gamma_0(asp_tmp) = gm_val_34(idx)*ones(m);
                end
            end
            
            %% 5. Drivetrain scaling:
            aspect_5 = 'Drivetrain';
            
            if(any(aspect_set == aspect_5))
                fprintf(opt_wrt, upper(aspect_5));
                
                fprintf('\n\tto be done later...\n\n');
            end
            
            %% Post processing:
            % Analysis of the residuals at each step:
            res_pp = [mean(res_1), mean(res_2),    res_3,    res_4];
            idx = ~isinf(res_pp);
            res_pp = res_pp(idx);
            
            asp_name = aspect_set;
            if(numel(aspect_set) ~= sum(idx))
                asp_name = asp_name(idx);
            end
            
            fprintf('Scale: %.1f kW = %.2f %% of Ref.\n', gamma_0('P')*[obj_ref.P_rated, ...
                                                                        100.0]);

            if(isempty(res_pp))
                fprintf('\tScaled version obtained by dimensional analysis scaling rules.\n');
            else
                [~, sorted_idx] = sort(res_pp);
                
                fprintf('\tOptimization residua:\n');

                for idx = sorted_idx
                    fprintf('\t%d. %s:\t%.5e\n', idx, asp_name{idx}, res_pp(idx));
                end
            end
            
            gamma_opt = gamma_0;
            
            res.stage = res_1;
            res.gear  = res_2;
            res.KJ    = res_3;
            res.KJ2   = res_4;
            res.SG    = res_12;
            res.KJg   = res_34;
            
            gamma_sep.stage = gm_val_1;
            gamma_sep.gear  = gm_val_2;
            gamma_sep.KJ    = gm_val_3;
            gamma_sep.KJ2   = gm_val_4;
            gamma_sep.SG    = gm_val_12;
            gamma_sep.KJg   = gm_val_34;
            
            gm0_cell = gamma_0.to_cell();
            obj_sca = scale_func(gm0_cell{:});
            
            for idx = 1:length(id)
                warning('on', id{idx});
            end
        end
        
        function [gamma, res, SH, f_n, mode_shape, k_mesh, gamma_asp, obj_array] = scaled_sweep(obj_ref, P_scale, varargin)
            %SCALED_SWEEP performs a sweep on the rated power parameter of
            % the object. Returns the scaling factors gamma
            %
            % see also SCALED_VERSION
            %
            
            try
                default = {'aspect_set', 'DA', ...
                    'n_scale', obj_ref.n_rotor};
                
                default = scaling_factor.process_varargin(default, varargin);
                
                fprintf('Reference %s drivetrain with rated power %.1f kW.\n', class(obj_ref), obj_ref.P_rated);
                
                n_P = numel(P_scale);
                n_fn = numel(obj_ref.resonances);
                
                gamma = zeros(length(obj_ref.gamma), n_P);
                res = struct;
                gamma_asp = struct;
                SH = zeros(numel(obj_ref.S_H), n_P);
                f_n = zeros(n_fn, n_P);
                mode_shape = zeros(n_fn, n_fn, n_P);
                k_mesh = zeros(numel(obj_ref.stage), n_P);
                
                gm_P = P_scale./obj_ref.P_rated;
                
                [~, MS_ref] = obj_ref.resonances;
                
                SH_ref = obj_ref.S_H;
                SS = [SH_ref, zeros(size(SH_ref))];
                
                plot_prop1 = {'lineStyle', '-' , 'lineWidth', 2.0, 'color', [346.6667e-3   536.0000e-3   690.6667e-3]};
                plot_prop2 = {'lineStyle', '--', 'lineWidth', 2.0, 'color', [915.2941e-3   281.5686e-3   287.8431e-3]};
                
                figure('units', 'centimeters', 'position', [5.0 5.0 34.0 12.0]);
                subplot(2, 6, 1:3)
                rectangle(obj_ref);
                xlim([0 6000])
                ylim([-1 1]*1500)
                title(sprintf('Reference: %.1f kW', obj_ref.P_rated));
                
                subplot(2, 6, 7:9)
                axis equal;
                xlim([0 6000])
                ylim([-1 1]*1500)
                
                for idx = 1:3
                    subplot(2, 6, idx + 3)
                    plot(1:14, MS_ref(:, idx)*90.0, plot_prop1{:});
                    
                    if(idx ~= 3)
                        subplot(2, 6, idx + 9)
                        plot(1:14, MS_ref(:, idx + 3)*90.0, plot_prop1{:})
                    else
                        subplot(2, 6, 12)
                        b = bar(SS);
                        ylim([1.0 2.25]);
                        yticks(1.0:0.25:2.25);
                        b(1).FaceColor = plot_prop1{end};
                        b(2).FaceColor = plot_prop2{end};
                    end
                end
                
                tick_x = ["R"  , ...
                    "c_1", "p_{11}", "p_{12}", "p_{13}", "s_1", ...
                    "c_2", "p_{21}", "p_{22}", "p_{23}", "s_2", ...
                    "W_3", "P_3", ...
                    "G"];
                idx_tick = 2:2:14;
                
                font_setting  = {'fontName', 'Times', 'fontSize', 12.0};
                latex_setting = {'fontName', 'Times', 'fontSize', 12.0, 'interpreter', 'LaTeX'};
                
                fig_axes = findobj(gcf, 'Type', 'Axes');
                
                set(get(fig_axes(1), 'ylabel'), 'string', '$S_H$'             , latex_setting{:});
                set(get(fig_axes(5), 'ylabel'), 'string', '$\theta_i$, [deg.]', latex_setting{:});
                set(get(fig_axes(6), 'ylabel'), 'string', '$\theta_i$, [deg.]', latex_setting{:});
                
                set(fig_axes(2:6), 'xlim'      , [1 14]);
                set(fig_axes(2:6), 'ylim'      , [-1 1]*100);
                set(fig_axes(2:6), 'xtick'     ,        idx_tick);
                set(fig_axes(2:6), 'xticklabel', tick_x(idx_tick));
                set(fig_axes     , font_setting{:});
                
                label_figure(fig_axes, font_setting);
                
                fig_name = @(i)(sprintf('@%s\\plots\\sweep_scale\\scale_%d_%d', class(obj_ref), i, n_P));
                
                file_name = sprintf('@%s\\plots\\sweep_scale\\scale_sweep.gif', class(obj_ref));
                
                k_P = mean(diff(P_scale));
                
                if(k_P > 0.0)
                    gamma_0 = scaling_factor(obj_ref.gamma.name, obj_ref.gamma.value*1.0e-3);
                else
                    gamma_0 = scaling_factor(obj_ref.gamma.name, obj_ref.gamma.value);
                end
                
                obj_array = cell(n_P + 1, 1);
                obj_array{1} = obj_ref;
                
                for idx = 1:n_P
                    [obj_sca, gamma_0, res_idx, gamma_idx] = obj_ref.scaled_version(P_scale(idx), ...
                        'gamma_0'   , gamma_0, ...
                        'n_scale'   , default.n_scale, ...
                        'aspect_set', default.aspect_set);
                    obj_array{idx + 1} = obj_sca;
                    
                    gamma(:, idx) = gamma_0.value;
                    res(idx).stage = res_idx.stage;
                    res(idx).gear  = res_idx.gear;
                    res(idx).KJ    = res_idx.KJ;
                    res(idx).KJ2   = res_idx.KJ2;
                    
                    gamma_asp(idx).stage = gamma_idx.stage;
                    gamma_asp(idx).gear  = gamma_idx.gear;
                    gamma_asp(idx).KJ    = gamma_idx.KJ;
                    gamma_asp(idx).KJ2   = gamma_idx.KJ2;
                    
                    SH(:, idx) = obj_sca.S_H;
                    k_mesh(:, idx) = [obj_sca.stage.k_mesh]';
                    [f_n(:, idx), mode_shape(:, :, idx)] = obj_sca.resonances;
                    
                    subplot(2, 6, 7:9);
                    cla;
                    rectangle(obj_sca);
                    xlim([0 6000]);
                    ylim([-1 1]*1500);
                    title(sprintf('Scale: %.1f kW = %.3f %% of Ref.', obj_sca.P_rated, gm_P(idx)*100.0));
                    SS(:, 2) = SH(:, idx);
                    
                    for jdx = 1:3
                        subplot(2, 6, jdx + 3);
                        cla;
                        hold on;
                        plot(1:14, MS_ref(:, jdx)*90.0, plot_prop1{:});
                        plot(1:14, mode_shape(:, jdx, idx)*90.0, plot_prop2{:});
                        
                        if(jdx ~= 3)
                            subplot(2, 6, jdx + 9);
                            cla;
                            hold on;
                            plot(1:14, MS_ref(:, jdx + 3)*90.0, plot_prop1{:})
                            plot(1:14, mode_shape(:, jdx + 3, idx)*90.0, plot_prop2{:});
                        else
                            subplot(2, 6, 12);
                            cla;
                            b = bar(SS);
                            hold on;
                            plot(0:10, 1.25*ones(1,11), 'k');
%                             yline(1.25, 'color', 'k');
                            ylim([1.0 2.25]);
                            yticks(1.0:0.25:2.25);
                            b(1).FaceColor = plot_prop1{end};
                            b(2).FaceColor = plot_prop2{end};
                        end
                    end
                    
                    set(get(fig_axes(1), 'ylabel'), 'string', '$S_H$, [-]'        , latex_setting{:});
                    set(fig_axes(1)  , 'xticklabel', ["s_1", "p_{1i}", "r_1", ...
                                                      "s_2", "p_{2i}", "r_2", ...
                                                      "P_3", "W_3"]);
                    
                    set(fig_axes, font_setting{:});
                    savefig(gcf, fig_name(idx));
                    print(fig_name(idx), '-dpng');
                    saveasGIF(file_name, idx);
                end
                
                close all;
            catch err
                save error_data_scaled_sweep;
                error(err.identifier, '%s', err.message);
            end
        end
        
        function [f, mode_shape] = scale_nth_resonance(obj, n, gamma_P, gamma_n, gamma, aspect)
            obj_sca = scale_aspect(obj, gamma_P, gamma_n, gamma, aspect);
            
            [f, mode_shape] = obj_sca.nth_resonance(n);
            
        end
        
        function [f_n, mode_shape] = scale_modal_analysis(obj, gamma_P, gamma_n, gamma, aspect)
            %SCALED_MODAL_ANALYSIS performs modal analysis on a scaled
            % Drivetrain object, returning only the N first resonances and
            % mode shapes.
            %
            % See also: MODAL_ANALYSIS.
            %
            
            obj_sca = scale_aspect(obj, gamma_P, gamma_n, gamma, aspect);
            
            [f_n, mode_shape] = obj_sca.resonances;
            
        end
        
        function [f_n, mode_shape] = scale_resonances(obj, gamma_val, gamma, N, normalize, aspect)
            %SCALED_RESONANCES returns the N first resonances and mode
            % shapes of a scaled Drivetrain object. The resonances can be
            % normalized or not.
            %
            
            lin = size(aspect, 1);
            
            for idx = 1:lin
                asp_tmp = aspect(idx, :);
                asp_tmp = asp_tmp(asp_tmp ~= '*');
                col = numel(asp_tmp);
                
                for jdx = 1:col
                    gamma(asp_tmp(jdx)) = gamma_val(idx, :);
                end
            end
            
            scale_func = str2func(class(obj));
            gm_cell = gamma.to_cell();
            obj_sca = scale_func(gm_cell{:});
            
            [f_n, mode_shape] = obj_sca.resonances('N', N, 'normalize', normalize);
            
        end
        
        function [SH, Sshaft] = scaled_safety_factors(obj_ref, gamma_P, gamma_n, gamma, aspect)
            %SCALED_SAFETY_FACTORS returns the safety factors of a scaled
            % Drivetrain object. Partial scaling is possible through the
            % argument opt_idx, which should be equal to [] for scaling all
            % parameters.
            %
            
            obj_sca = obj_ref.scale_aspect(gamma_P, gamma_n, gamma, aspect);
            
            SH     = obj_sca.S_H;
            Sshaft = obj_sca.S_shaft;

        end
        
    end
    
    %% Set methods:
    methods
        function obj = set.dynamic_model(obj, val)
            obj.dynamic_model = val;
        end
        
    end
    
    %% Get methods:
    methods
        function val = get.n_out(obj)
            val = zeros(obj.N_stage, 1);
            
            for idx = 1:obj.N_stage
                val(idx) = obj.n_rotor*prod([obj.stage(1:idx).u]);
            end
        end
        
        function val = get.T_out(obj)
            val = (obj.P_rated*1.0e3)./(obj.n_out*pi/30.0);
        end
        
        function val = get.u(obj)
            val = prod([obj.stage(1:end).u]);
        end
        
        function val = get.S_H(obj)
            val = obj.S_H_val;
            
            if(isrow(val))
                val = val';
            end
        end
        
        function val = get.S_F(obj)
            val = obj.S_F_val;
            
            if(isrow(val))
                val = val';
            end
        end
        
        function val = get.S_shaft(obj)
            val = obj.S_shaft_val;
        end
        
        function val = get.S_G(obj)
            val = [obj.S_H;
                obj.S_F];
            val(isnan(val)) = [];
        end
            
    end
    
    %% Misc
    methods(Static)
		function val = read_if2(file_name)
			%READ_IF2 converts .if2 files to .mat files.
				
			file_ID = fopen(file_name, 'r');

			file_type_descriptor = fscanf(file_ID, '%s', 1);
			if(~strcmp(file_type_descriptor, '$SIMPACK_Input_Function_Set$'))
				error('file type wrong.')
			end

			release = fscanf(file_ID, '%f ! Release \n', 1);
			file_format = fscanf(file_ID, '%f ! Format: %*s = %*s\n');

			name = fscanf(file_ID, '%s ! Input Function 1\n', 1);
			interp_method  = fscanf(file_ID, ' %f !            Interpolation Method: %*s = %*s \n', 1);
			extrapol_param = fscanf(file_ID, '  %f !            Extrapolation: %*s = %*s ; %*s\n', [1 2]);
			scaling_factor = fscanf(file_ID, '  %f \t %f !            Unit Factors %*s = %*s', [1 2]);
			unit_fac1 = fscanf(file_ID, '%s', 1);
			unit_fac2 = fscanf(file_ID, ']  | [%s] \n', 1);
			unit_typ1 = fscanf(file_ID, '%s !            Unit Types \n', 1);
			unit_typ2 = fscanf(file_ID, '%s !            Unit Types \n', 1);

			data = fscanf(file_ID, '%f', [2, Inf])';

			fclose(file_ID);
			
			val.release = release;
			val.file_format = file_format;
			val.name = name;
			val.interp_method = interp_method;
			val.extrapol_param = extrapol_param;
			val.scaling_factor = scaling_factor;
			val.unit_type = [string(unit_typ1), string(unit_typ2)];
			val.time = data(:,1);
			val.data = data(:,2);

        end
        
        function rewrite_subvar(old_file)
            %REWRITE_SUBVAR processes a .subvar file making it easier to
            % modify later on using the method update_subvar() that should 
            % be implemented by the sub-classes.
            %
            % see also update_subvar
            %
            
            old_ID = fopen(old_file, 'r');
            new_ID = fopen('new.subvar', 'w');

            while ~feof(old_ID)
                old_line = fgetl(old_ID);
                fprintf(new_ID, 'fprintf(new_ID, ''%s\\n'');\n', old_line);
            end
            fclose(old_ID);
            fclose(new_ID);
        end
    end
end
