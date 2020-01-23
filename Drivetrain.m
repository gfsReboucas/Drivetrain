classdef Drivetrain
    %DRIVETRAIN This class implements SOME procedures for the dynamic
    % analysis and scaling of drivetrains. The safety factor for surface 
    % durability (pitting) is calculated according to ISO 6336 [1, 2]. The
    % NREL 5MW reference gearbox proposed by Nejad et. al. [3] is used as
    % the default drivetrain model, but other models should be implemented
    % in the future.
    %
    % References:
    % [1] ISO 6336-1:2006 Calculation of load capacity of spur and helical
    % gears -- Part 1: Basic principles, introduction and general influence
    % factors
    % [2] ISO 6336-2:2006 Calculation of load capacity of spur and helical
    % gears -- Part 2: Calculation of surface durability (pitting)
    % [3] Nejad, A. R., Guo, Y., Gao, Z., Moan, T. (2016). Development of a
    % 5 MW reference gearbox for offshore wind turbines. Wind Energy. 
    % https://doi.org/10.1002/we.1884
    % [4] Jonkman, J., Butterfield, S., Musial, W., & Scott, G. Definition
    % of a 5-MW Reference Wind Turbine for Offshore System Development. 
    % doi:10.2172/947422.
    % [5] IEC 61400-4:2012 Wind Turbines -- Part 4: Design Requirements for
    % wind turbine gearboxes.
    % [6] Anaya-Lara, O., Tande, J.O., Uhlen, K., Merz, K. and Nejad, A.R.
    % (2019). Modelling and Analysis of Drivetrains in Offshore Wind
    % Turbines. In Offshore Wind Energy Technology (eds O. Anaya-Lara, J.O.
    % Tande, K. Uhlen and K. Merz). doi:10.1002/9781119097808.ch3
    %
    % written by:
    % Geraldo Rebou�as
    % - Geraldo.Reboucas@ntnu.no OR
    % - gfs.reboucas@gmail.com
    %
    % Postdoctoral Fellow at:
    %    Norwegian University of Science and Technology, NTNU
    %    Department of Marine Technology, IMT
    %    Marine System Dynamics and Vibration Lab, MD Lab
    %    https://www.ntnu.edu/imt/lab/md-lab
    %
    
    properties(SetAccess = private)
        stage       (1, :) Gear_Set;                                                            % [-],      gearbox stages
        P_rated     (1, :)          {mustBeNumeric, mustBeFinite, mustBePositive} = 5.0e3;      % [kW],     Rated power
        n_rotor     (1, :)          {mustBeNumeric, mustBeFinite, mustBePositive} = 12.1;       % [1/min.], Rated rotor speed
        main_shaft  (1, :) Shaft                                                  = Shaft;      % [-],      Input Shaft
        m_Rotor     (1, :)          {mustBeNumeric, mustBeFinite, mustBePositive} = 110.0e3;    % [kg],     Rotor mass according to [3]
        J_Rotor     (1, :)          {mustBeNumeric, mustBeFinite, mustBePositive} = 57231535.0; % [kg-m^2], Rotor mass moment of inertia according to [6]
        m_Gen       (1, :)          {mustBeNumeric, mustBeFinite, mustBePositive} = 1900.0;     % [kg],     Generator mass according to [4]
        J_Gen       (1, :)          {mustBeNumeric, mustBeFinite, mustBePositive} = 534.116;    % [kg-m^2], Generator mass moment of inertia [4]
        N_stg       (1, 1)          {mustBeNumeric, mustBeFinite, mustBePositive} = 3;          % [-],      Number of stages
    end
    
    properties(SetAccess = private)
        S_shaft_val (1, 4)          {mustBeNumeric, mustBeFinite, mustBePositive} = 1.0;    % [-],      Safey factor for the shafts
        S_H_val     (1, 6)          {mustBeNumeric, mustBeFinite, mustBePositive} = 1.25;   % [-],      Safety factor for surface durability (against pitting)
    end
    
    properties
        dynamic_model (1, :) string {mustBeMember(dynamic_model, ["Thomson_ToV", ...
                                                                  "Kahraman_1994", ...
                                                                  "Lin_Parker_1999"])} = "Thomson_ToV"; % which dynamic model should be used to perform modal analysis on the Drivetrain.
    end
    
    properties(Dependent)
        T_out;     % [N-m],    Output torque for each stage
        n_out;     % [1/min.], Output speed  for each stage
        u;       % [-],      Cumulative gear ratio
        S_H;     % [-],      Safety factor for surface durability (against pitting)
        S_shaft; % [-],      Safey factor for the shafts
    end
    
    methods
        function obj = Drivetrain(N_st, stage, P_r, n_r, inp_shaft, m_R, J_R, m_G, J_G)
            if(nargin == 0)
                N_st = 3;
                for idx = 1:N_st
                    stage(idx) = Gear_Set.NREL_5MW(idx);
                end
                
                P_r = 5.0e3; % [kW], Rated power
                n_r = 12.1; % [1/min.], Input speed
                inp_shaft = Shaft;
                
                m_R = 110.0e3;      J_R = 57231535.0;
                m_G = 1900.0;       J_G = 534.116;
                
            end
            
            obj.stage = stage;
            obj.P_rated = P_r;
            obj.n_rotor = n_r;
            obj.main_shaft = inp_shaft;
            
            obj.m_Rotor = m_R;          obj.J_Rotor = J_R;
            obj.m_Gen = m_G;            obj.J_Gen = J_G;
            
            obj.N_stg = N_st;
            
            [obj.S_H_val, obj.S_shaft_val] = obj.safety_factors;
            
        end
        
        function tab = disp(obj)
            %DISP display some properties of a Drivetrain object
            % description, symbol, unit, value
            
            tab_str = {"Rated power",                                  "P",       "kW",     "-+-",                obj.P_rated,          "-+-";            % 1
                       "Output Speed (Sun/Pinion)",                    "n_out",   "1/min.", obj.n_out(1),         obj.n_out(2),         obj.n_out(3);
                       "Output Torque (Sun/Pinion)",                   "T_out",   "N-m",    obj.T_out(1),         obj.T_out(2),         obj.T_out(3);
                       "Minimum safety factor against pitting",        "S_Hmin",  "-",      1.25,                 1.25,                 1.25;
                       "Safety factor against pitting (Sun/Pinion)",   "S_H1",    "-",      obj.S_H(1),           obj.S_H(3),           obj.S_H(5);
                       "Safety factor against pitting (Planet/Wheel)", "S_H2",    "-",      obj.S_H(2),           obj.S_H(4),           obj.S_H(6);
                       "Safety factor (Shaft)",                        "S",        "-",     obj.S_shaft(2),       obj.S_shaft(3),       obj.S_shaft(4);
                       "Type",                                         "-",       "-",      obj.stage(:).configuration;
                       "Gear ratio",                                   "u",       "-",      obj.stage(:).u;
                       "Number of planets",                            "p",       "-",      obj.stage(:).N_p;
                       "Normal module",                                "m_n",     "mm",     obj.stage(:).m_n;
                       "Normal pressure angle",                        "alpha_n", "deg.",   obj.stage(:).alpha_n;
                       "Helix angle",                                  "beta",    "deg.",   obj.stage(:).beta;
                       "Face width",                                   "b",       "mm",     obj.stage(:).b;
                       "Center distance",                              "a_w",     "mm",     obj.stage(:).a_w;
                       "Number of teeth (Sun/Pinion)",                 "z_1",     "-",      obj.stage(1).z(1),    obj.stage(2).z(1),    obj.stage(3).z(1);
                       "Number of teeth (Planet/Wheel)",               "z_2",     "-",      obj.stage(1).z(2),    obj.stage(2).z(2),    obj.stage(3).z(2);
                       "Number of teeth (Ring)",                       "z_3",     "-",      obj.stage(1).z(3),    obj.stage(2).z(3),    "-*-";
                       "Profile shift coefficient (Sun/Pinion)",       "x_1",     "-",      obj.stage(1).x(1),    obj.stage(2).x(1),    obj.stage(3).x(1);
                       "Profile shift coefficient (Planet/Wheel)",     "x_2",     "-",      obj.stage(1).x(2),    obj.stage(2).x(2),    obj.stage(3).x(2);
                       "Profile shift coefficient (Ring)",             "x_3",     "-",      obj.stage(1).x(3),    obj.stage(2).x(3),    "-*-";
                       "Reference diameter (Sun/Pinion)",              "d_1",     "mm",     obj.stage(1).d(1),    obj.stage(2).d(1),    obj.stage(3).d(1);
                       "Reference diameter (Planet/Wheel)",            "d_2",     "mm",     obj.stage(1).d(2),    obj.stage(2).d(2),    obj.stage(3).d(2);
                       "Reference diameter (Ring)",                    "d_3",     "mm",     obj.stage(1).d(3),    obj.stage(2).d(3),    "-*-";
                       "Tip diameter (Sun/Pinion)",                    "d_a1",    "mm",     obj.stage(1).d_a(1),  obj.stage(2).d_a(1),  obj.stage(3).d_a(1);
                       "Tip diameter (Planet/Wheel)",                  "d_a2",    "mm",     obj.stage(1).d_a(2),  obj.stage(2).d_a(2),  obj.stage(3).d_a(2);
                       "Tip diameter (Ring)",                          "d_a3",    "mm",     obj.stage(1).d_a(3),  obj.stage(2).d_a(3),  "-*-";
                       "Root diameter (Sun/Pinion)",                   "d_f1",    "mm",     obj.stage(1).d_f(1),  obj.stage(2).d_f(1),  obj.stage(3).d_f(1);
                       "Root diameter (Planet/Wheel)",                 "d_f2",    "mm",     obj.stage(1).d_f(2),  obj.stage(2).d_f(2),  obj.stage(3).d_f(2);
                       "Root diameter (Ring)",                         "d_f3",    "mm",     obj.stage(1).d_f(3),  obj.stage(2).d_f(3),  "-*-";
                       "Mass (Sun/Pinion)",                            "m_1",     "kg",     obj.stage(1).mass(1), obj.stage(2).mass(1), obj.stage(3).mass(1);
                       "Mass (Planet/Wheel)",                          "m_2",     "kg",     obj.stage(1).mass(2), obj.stage(2).mass(2), obj.stage(3).mass(2);
                       "Mass (Ring)",                                  "m_3",     "kg",     obj.stage(1).mass(3), obj.stage(2).mass(3), "-*-";
                       "Mass moment of inertia (Sun/Pinion)",          "J_xx1",   "kg-m^2", obj.stage(1).J_x(1),  obj.stage(2).J_x(1),  obj.stage(3).J_x(1);
                       "Mass moment of inertia (Planet/Wheel)",        "J_xx2",   "kg-m^2", obj.stage(1).J_x(2),  obj.stage(2).J_x(2),  obj.stage(3).J_x(2);
                       "Mass moment of inertia (Ring)",                "J_xx3",   "kg-m^2", obj.stage(1).J_x(3),  obj.stage(2).J_x(3),  "-*-";
                       "Mass moment of inertia (Sun/Pinion)",          "J_yy1",   "kg-m^2", obj.stage(1).J_y(1),  obj.stage(2).J_y(1),  obj.stage(3).J_y(1);
                       "Mass moment of inertia (Planet/Wheel)",        "J_yy2",   "kg-m^2", obj.stage(1).J_y(2),  obj.stage(2).J_y(2),  obj.stage(3).J_y(2);
                       "Mass moment of inertia (Ring)",                "J_yy3",   "kg-m^2", obj.stage(1).J_y(3),  obj.stage(2).J_y(3),  "-*-";
                       "Mass moment of inertia (Sun/Pinion)",          "J_zz1",   "kg-m^2", obj.stage(1).J_z(1),  obj.stage(2).J_z(1),  obj.stage(3).J_z(1);
                       "Mass moment of inertia (Planet/Wheel)",        "J_zz2",   "kg-m^2", obj.stage(1).J_z(2),  obj.stage(2).J_z(2),  obj.stage(3).J_z(2);
                       "Mass moment of inertia (Ring)",                "J_zz3",   "kg-m^2", obj.stage(1).J_z(3),  obj.stage(2).J_z(3),  "-*-";
                       };

            Parameter = tab_str(:,1);
            Symbol    = tab_str(:,2);
            Unit      = tab_str(:,3);
            Stage_01  = tab_str(:,4);
            Stage_02  = tab_str(:,5);
            Stage_03  = tab_str(:,6);

            tab = table(           Parameter,   Symbol,   Stage_01,   Stage_02,   Stage_03,   Unit, ...
                'variableNames', ["Parameter", "Symbol", "Stage_01", "Stage_02", "Stage_03", "Unit"]);
            
            if(nargout == 0)
                fprintf("Gear box:\n");
                disp(tab);
                fprintf("Main shaft:\n");
                disp(obj.main_shaft);
                clear tab;
            end

        end
        
        function [tab, tab_str] = comp_stage(ref, sca, idx)
            tab_str = {"Rated power",                                  "P",       "kW",     ref.P_rated,                sca.P_rated,                sca.P_rated               /ref.P_rated;                % 1
                       "Output Speed (Sun/Pinion)",                    "n_out",   "1/min.", ref.n_out(idx),             sca.n_out(idx),             sca.n_out(idx)            /ref.n_out(idx);             % 2
                       "Output Torque (Sun/Pinion)",                   "T_out",   "N-m",    ref.T_out(idx),             sca.T_out(idx),             sca.T_out(idx)            /ref.T_out(idx);             % 3
                       "Safety factor against pitting (Sun/Pinion)",   "S_H1",    "-",      ref.S_H(2*idx - 1),         sca.S_H(2*idx - 1),         sca.S_H(2*idx - 1)        /ref.S_H(2*idx - 1);         % 4
                       "Safety factor against pitting (Planet/Wheel)", "S_H2",    "-",      ref.S_H(2*idx),             sca.S_H(2*idx),             sca.S_H(2*idx)            /ref.S_H(2*idx);             % 5
                       "Gear ratio",                                   "u",       "-",      ref.stage(idx).u,           sca.stage(idx).u,           sca.stage(idx).u          /ref.stage(idx).u;           % 6
                       "Normal module",                                "m_n",     "mm",     ref.stage(idx).m_n,         sca.stage(idx).m_n,         sca.stage(idx).m_n        /ref.stage(idx).m_n;         % 7
                       "Face width",                                   "b",       "mm",     ref.stage(idx).b,           sca.stage(idx).b,           sca.stage(idx).b          /ref.stage(idx).b;           % 8
                       "Center distance",                              "a_w",     "mm",     ref.stage(idx).a_w,         sca.stage(idx).a_w,         sca.stage(idx).a_w        /ref.stage(idx).a_w;         % 9
                       "Reference diameter (Sun/Pinion)",              "d_1",     "mm",     ref.stage(idx).d(1),        sca.stage(idx).d(1),        sca.stage(idx).d(1)       /ref.stage(idx).d(1);        % 10
                       "Mass (Sun/Pinion)",                            "m_1",     "kg",     ref.stage(idx).mass(1),     sca.stage(idx).mass(1),     sca.stage(idx).mass(1)    /ref.stage(idx).mass(1);     % 11
                       "Mass moment of inertia (Sun/Pinion)",          "J_xx1",   "kg-m^2", ref.stage(idx).J_x(1),      sca.stage(idx).J_x(1),      sca.stage(idx).J_x(1)     /ref.stage(idx).J_x(1);      % 12
                       "Shaft / Diameter",                             "d",       "mm",     ref.stage(idx).out_shaft.d, sca.stage(idx).out_shaft.d, sca.stage(idx).out_shaft.d/ref.stage(idx).out_shaft.d; % 13
                       "Shaft / Length",                               "L",       "mm",     ref.stage(idx).out_shaft.L, sca.stage(idx).out_shaft.L, sca.stage(idx).out_shaft.L/ref.stage(idx).out_shaft.L; % 14
                        };

            Reference = tab_str(:, 4);
            Scale     = tab_str(:, 5);
            Ratio     = tab_str(:, 6);
            
            tab = table(Scale, Reference, Ratio, ...
                    'VariableNames', ["Scale", "Reference", "Ratio"]);
        end
        
        function tab = disp_comp(ref, sca)
             Stage_01           = comp_stage(ref, sca, 1);
             Stage_02           = comp_stage(ref, sca, 2);
            [Stage_03, tab_str] = comp_stage(ref, sca, 3);
            
            Parameter = tab_str(:,1);
            Symbol    = tab_str(:,2);
            Unit      = tab_str(:,3);
            
            tab = table(           Parameter,   Symbol,   Stage_01,   Stage_02,   Stage_03,   Unit, ...
                'variableNames', ["Parameter", "Symbol", "Stage_01", "Stage_02", "Stage_03", "Unit"]);
            disp(tab);
            
            if(nargout == 0)
                clear tab;
            end

        end
        
        function plot(obj)
            hold on;
            axis equal;
            box on;
            for idx = 1:obj.N_stg
                subplot(1, obj.N_stg, idx)
                obj.stage(idx).plot;
            end
        end
        
        function rectangle(obj)
            
            % LINSPECER: Plot lots of lines with very distinguishable and 
            % aesthetically pleasing colors. It can be dowloaded from
            % MATLAB's File Exchange on:
            % https://se.mathworks.com/matlabcentral/fileexchange/42673-beautiful-and-distinguishable-line-colors-+-colormap
            color = linspecer(6, "qualitative");
            
            hold on;
            C_s = [obj.main_shaft.L/2.0 0.0]';
            rectangle(obj.main_shaft, C_s, color(6, :), "edgeColor", "k", "lineStyle", "-" , "faceColor", color(6, :));
            
            C = [obj.main_shaft.L;
                 obj.stage(1).carrier.b + obj.stage(1).out_shaft.L;
                 obj.stage(2).carrier.b + obj.stage(2).out_shaft.L];

            for idx = 1:obj.N_stg
                rectangle(obj.stage(idx), [sum(C(1:idx)) 0.0]');
            end
            hold off;
        end
        
        function plot_comp(DT1, DT2)
            if(DT1.N_stg ~= DT2.N_stg)
                error("Both DT's should have the same number of stages.");
            end
            
            hold on;
            axis equal;
            box on;
            for idx = 1:DT1.N_stg
                subplot(2, DT1.N_stg, idx)
                DT1.stage(idx).plot;
                title(sprintf("Stage %d", idx));
                subplot(2, DT2.N_stg, idx + DT2.N_stg)
                DT2.stage(idx).plot;
            end

        end
        
    end
    
    methods
        %% Calculations:
        function obj_sca = scale_Drivetrain(obj_ref, gamma_P, gamma_n, gamma)
            %SCALE_DRIVETRAIN returns a Drivetrain object scaled by the
            % factors gamma_P and gamma_n for its rated power and rotor
            % speed and gamma for normal module, face width, shaft
            % dimensions (diameter and length), mass and mass moment of
            % inertia of rotor and generator.
            %
            % current limitations/constraints:
            % - shaft diameter is proportional to torque;
            % - dimensions of the main shaft are proportional to the
            % output shaft in stage 03;
            %
            
            if((length(gamma) < 18) || (~isa(gamma, "containers.Map")))
                error("gamma must be a container with 18 elements.");
            end
            
            %          Normal module, Face width, Length, Diameter (output shaft)
            key_set = ["m_n1"       , "b_1"     , "L_1" , "d_1", ... % Stage 01
                       "m_n2"       , "b_2"     , "L_2" , "d_2", ... % Stage 02
                       "m_n3"       , "b_3"     , "L_3" , "d_3", ... % Stage 03
               ... %   Mass         , M. M. Iner, Mass  , M. M. Iner,
                       "m_R"        , "J_R"     , "m_G" , "J_G", ... % Rotor and Generator
               ... %   Length       , Diameter
                       "L_s"        , "d_s"];        %#ok<*CLARRSTR> % Main shaft
            
            gamma_mR  = gamma("m_R");
            gamma_JR  = gamma("J_R");
            gamma_mG  = gamma("m_G");
            gamma_JG  = gamma("J_G");
%             gamma_Ls  = gamma("L_s");
%             gamma_ds  = gamma("d_s");
            
            gamma_Torque = gamma_P/gamma_n; % Applied torque
            
            % Main shaft:
%             gamma_Ls = gamma("L_3"); % length
            gamma_ds = nthroot(gamma_Torque, 3.0);
            
            % Shaft diameters are scaled in the same way
            gamma("d_1") = gamma_ds;
            gamma("d_2") = gamma_ds;
            gamma("d_3") = gamma_ds;
            gamma("d_s") = gamma_ds;
            
            stage_sca = [Gear_Set, Gear_Set, Gear_Set];
            
            idx_mn = find(contains(key_set, "m_n"));
            idx_bb = find(contains(key_set, "b"));
            idx_LL = find(contains(key_set, "L"));
            idx_dd = find(contains(key_set, "d"));
            
            gamma_mn = cell2mat(values(gamma, {key_set{idx_mn}}))'; %#ok<*FNDSB>
            gamma_bb = cell2mat(values(gamma, {key_set{idx_bb}}))';
            gamma_LL = cell2mat(values(gamma, {key_set{idx_LL}}))';
            gamma_dd = cell2mat(values(gamma, {key_set{idx_dd}}))';
            
            for idx = 1:obj_ref.N_stg
                % Scaled stage idx:
                gamma_stage = [gamma_mn(idx) gamma_bb(idx)]';
                % Scaled output shaft idx:
                gamma_shaft = [gamma_dd(idx) gamma_LL(idx)]';
                stage_sca(idx) = obj_ref.stage(idx).scaled_Gear_Set(gamma_stage, gamma_shaft);
            end
            
            idx = idx + 1;
            
            main_shaft_sca = Shaft(obj_ref.main_shaft.d*gamma_dd(idx), ...
                                   obj_ref.main_shaft.L*gamma_LL(idx));
            
            obj_sca = Drivetrain(obj_ref.N_stg, ...
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
        
        function obj_sca = scale_aspect(obj_ref, gamma_P, gamma_n, gamma, aspect)
            %SCALE_ASPECT returns a Drivetrain object scaled by the factors
            % gamma_P and gamma_n for its rated power and rotor speed. The
            % argument gamma can be used to scale the Drivetrain
            % considering various aspects.
            %
            
            %          Normal module, Face width, Length, Diameter (output shaft)
            key_set = ["m_n1"       , "b_1"     , "L_1" , "d_1", ... % Stage 01
                       "m_n2"       , "b_2"     , "L_2" , "d_2", ... % Stage 02
                       "m_n3"       , "b_3"     , "L_3" , "d_3", ... % Stage 03
               ... %   Mass         , M. M. Iner, Mass  , M. M. Iner,
                       "m_R"        , "J_R"     , "m_G" , "J_G", ... % Rotor and Generator
               ... %   Length       , Diameter
                       "L_s"        , "d_s"];                        % Main shaft
            
            gamma_full = containers.Map(key_set, ones(18, 1));
            
            switch(aspect)
                case "Drivetrain"
                    gamma_full = containers.Map(key_set, gamma);
                    
                case "by_stage"
                    % scaled by stage (i.e. one scale factor for both the
                    % normal module and face width per stage).
                    
                    if(numel(gamma) < 3)
                        error("gamma must have 3 elements.");
                    end
                    
                    sub_key ={"m_n1", "b_1", ...
                              "m_n2", "b_2", ...
                              "m_n3", "b_3"};
                          
                    for idx = 1:3
                        gamma_full(sub_key{2*idx - 1}) = gamma(idx, :);
                        gamma_full(sub_key{2*idx})     = gamma(idx, :);
                    end

                case "gear"
                    % scales only the gear dimensions (normal module and 
                    % face width).
                    
                    if(numel(gamma) < 6)
                        error("gamma must have 6 elements.");
                    end

                    sub_key ={"m_n1", "b_1", ...
                              "m_n2", "b_2", ...
                              "m_n3", "b_3"};
                          
                    for idx = 1:6
                        gamma_full(sub_key{idx}) = gamma(idx, :);
                    end

                case "stiff_mass_mom_inertia"
                    % scales only parameters related to the shaft's 
                    % stiffness (length) and mass moment of inertia of 
                    % rotor and generator.

                    if(numel(gamma) < 2)
                        error("gamma must have 2 elements.");
                    end
                    
                    sub_key = {"L_1", "L_2", "L_3", ...
                               "J_R", "J_G", ...
                               "L_s"};
                          
                    for idx = 1:3
                        gamma_full(sub_key{idx}) = gamma(1, :);
                    end
                    
                    for idx = 4:5
                        gamma_full(sub_key{idx}) = gamma(2, :);
                    end
                    
                    gamma_full(sub_key{6}) = gamma(1, :);
                    
                case "stiff_mass_mom_inertia_detail"
                    % scales only parameters related to the shaft's 
                    % stiffness (length) and mass moment of inertia of 
                    % rotor and generator.

                    if(numel(gamma) < 6)
                        error("gamma must have 6 elements.");
                    end
                    
                    sub_key = {"L_1", "L_2", "L_3", ...
                               "J_R", "J_G", ...
                               "L_s"};
                          
                    for idx = 1:6
                        gamma_full(sub_key{idx}) = gamma(idx, :);
                    end
                    
                case "shaft_stiffness"
                    % scales only parameters related to the shaft's 
                    % stiffness (length).
                    
                    if(numel(gamma) < 3)
                        error("gamma must have 3 elements.");
                    end
                    
                    idx_L = find(contains(key_set, "L"));
                    
                    for idx = 1:3
                        gamma_full(key_set{idx_L(idx)}) = gamma(idx, :);
                    end

                case "stiffness"
                    % scales only stiffness-related parameters. The mesh
                    % stiffness is assumed to be proportinal to the gear's 
                    % face width.
                    
                    if(numel(gamma) < 6)
                        error("gamma must have 6 elements.");
                    end

                    gamma_full(2,  :) = gamma(1, :); % face width, stage 01 
                    gamma_full(3,  :) = gamma(2, :); % shaft length
                    gamma_full(6,  :) = gamma(3, :); % face width, stage 02
                    gamma_full(7,  :) = gamma(4, :); % shaft length
                    gamma_full(10, :) = gamma(5, :); % face width, stage 03
                    gamma_full(11, :) = gamma(6, :); % shaft length

                otherwise
                    error("Option [%s] is NOT valid.", aspect);
            end
            
            obj_sca = obj_ref.scale_Drivetrain(gamma_P, gamma_n, gamma_full);
            
        end
           
        function [obj_sca, gamma, res, gamma_sep] = scaled_version(obj_ref, P_scale, n_R_scale, normalize_freq, N_freq, varargin)
            %SCALED_VERSION returns an scaled Drivetrain object together
            % with the scaling factor for its parameters and residual error
            % from the scaling optimization process. Different scaling can
            % be obtained by using:
            % - model_approach: Different model approaches;
            % - normalize_freq: Normalizing or not the target resonances;
            % - N_freq: choosing how many resonances should be considered
            % during the scaling optimization process;
            % Additionaly, one can define an initial value for gamma to be
            % used on the optimization process
            
            if(~isempty(varargin))
                gamma_prev = varargin{1};
            else
                gamma_prev = ones(18, 1)*0.5;
            end
            
            % Input scaling factors:
            gamma_P = P_scale/obj_ref.P_rated;
            gamma_n = n_R_scale/obj_ref.n_rotor;
            
            %% 1. Stage scaling:
            aspect_01 = "by_stage";
            
            S_H_ref = obj_ref.S_H;
            
            % function for aspect_01:
            fun_asp  = @(x)(S_H_ref - obj_ref.scaled_safety_factors(gamma_P, gamma_n, x, aspect_01));
            fun_min = @(x)(norm(fun_asp(x))^2);
            
            gamma_min = ones(3, 1)*1.0e-6;
            gamma_Max = ones(3, 1);
%             gamma_0 = gamma_Max*nthroot(gamma_P, 3.0); % for faster convergence
            gamma_0 = gamma_prev(1:2:5);
            
            fun_ineq = @(x)(1.25 - obj_ref.scaled_safety_factors(gamma_P, gamma_n, x, aspect_01));

            constraint_fun = @(x)deal(fun_ineq(x), fun_asp(x)); % inequalities, equalities
            
            opt_solver = optimoptions("fmincon", "display", "notify");
            
            id_1 = "prog:input";
            id_2 = "MATLAB:nearlySingularMatrix";
            warning("off", id_1);
            warning("off", id_2);
            [gamma_01, res_01] = fmincon(fun_min, gamma_0, [], [], [], [], gamma_min, gamma_Max, constraint_fun, opt_solver);
            
            %% 2. Gears:
            aspect_02 = "gear";
            
            % function for aspect_02:
            fun_asp  = @(x)(S_H_ref - obj_ref.scaled_safety_factors(gamma_P, gamma_n, x, aspect_02));
            fun_min = @(x)(norm(fun_asp(x))^2);
            
            gamma_min = ones(6, 1)*1.0e-6;
            gamma_Max = ones(6, 1);
            
            gamma_0 = gamma_01([1 1 2 2 3 3]);
%             gamma_0 = gamma_prev(1:6);
            
            fun_ineq = @(x)(1.25 - obj_ref.scaled_safety_factors(gamma_P, gamma_n, x, aspect_02));

            constraint_fun = @(x)deal(fun_ineq(x), fun_asp(x)); % inequalities, equalities
            
            [gamma_02, res_02, ~] = fmincon(fun_min, gamma_0, [], [], [], [], gamma_min, gamma_Max, constraint_fun, opt_solver);
            
            obj_sca_02 = obj_ref.scale_aspect(gamma_P, gamma_n, gamma_02, aspect_02);
            
            %% 3. Shaft stiffness and Mass moment of inertia of rotor and generator:
            aspect_03 = "stiff_mass_mom_inertia";
            
            f_n_ref = obj_ref.resonances(N_freq, normalize_freq);
            
            % function for aspect_03: gamma_P and gamma_n are set to 1.0
            % because the Drivetrain is already scaled for these
            % parameters.
            fun_asp  = @(x)(1.0 - obj_sca_02.scale_resonances(N_freq, normalize_freq, 1.0, 1.0, x, aspect_03)./f_n_ref);
            fun_min = @(x)(norm(fun_asp(x))^2);

            gamma_min = ones(2, 1)*1.0e-6;
            gamma_Max = ones(2, 1);
            
%             gamma_0 = [ones(3, 1)*nthroot(gamma_P, 3.0);
%                        ones(2, 1)*  power(gamma_P, 1.6)];
%             gamma_0 = gamma_Max*0.5;
            gamma_0 = gamma_prev([3 14]);
            
            constraint_fun = @(x)deal([], fun_asp(x)); % inequalities, equalities
            
            [gamma_03, res_03, ~] = fmincon(fun_min, gamma_0, [], [], [], [], gamma_min, gamma_Max, constraint_fun, opt_solver);

            obj_sca = obj_sca_02.scale_aspect(1.0, 1.0, gamma_03, aspect_03);
            
            %% 4. Detailed version of 3:
%             aspect_04 = "stiff_mass_mom_inertia_detail";
%             
%             % function for aspect_03: gamma_P and gamma_n are set to 1.0
%             % because the Drivetrain is already scaled for these
%             % parameters.
%             fun_asp  = @(x)(1.0 - obj_sca_02.scale_resonances(N_freq, normalize_freq, 1.0, 1.0, x, aspect_04)./f_n_ref);
%             fun_min = @(x)(norm(fun_asp(x))^2);
% 
%             gamma_min = ones(6, 1)*1.0e-6;
%             gamma_Max = ones(6, 1);
%             
% %             gamma_0 = [ones(3, 1)*nthroot(gamma_P, 3.0);
% %                        ones(2, 1)*  power(gamma_P, 1.6)];
% %             gamma_0 = gamma_Max*0.5;
%             gamma_0 = gamma_03([1 1 1 2 2 1]);
%             
%             constraint_fun = @(x)deal([], fun_asp(x)); % inequalities, equalities
%             
%             [gamma_04, res_04, ~] = fmincon(fun_min, gamma_0, [], [], [], [], gamma_min, gamma_Max, constraint_fun, opt_solver);
%             
%             obj_sca = obj_sca_02.scale_aspect(1.0, 1.0, gamma_04, aspect_04);
            warning("on", id_2);
            warning("on", id_1);
            
            %% Post processing:
            % Analysis of the residuals at each step:
            res_pp = [res_01, res_02, res_03];
            aspect = [aspect_01, aspect_02, aspect_03];
            
            [~, sorted_idx] = sort(res_pp);
            
            fprintf("Scale: %.1f kW = %.2f %% of Ref.\n", obj_sca.P_rated, gamma_P*100.0);
            for idx = sorted_idx
                fprintf("%d. %s:\t%.5e\n", idx, aspect(idx), res_pp(idx));
            end
            
            gamma = [gamma_02(1), gamma_02(2), gamma_03(1), 1.0        , ...
                     gamma_02(3), gamma_02(4), gamma_03(1), 1.0        , ...
                     gamma_02(5), gamma_02(6), gamma_03(1), 1.0        , ...
                     1.0        , gamma_03(2), 1.0        , gamma_03(2), ...
                     gamma_03(1), 1.0]';
                 
            res.stage = res_01;
            res.gear = res_02;
            res.KJ = res_03;
%             res.KJ2 = res_04;
            
            gamma_sep.stage = gamma_01;
            gamma_sep.gear = gamma_02;
            gamma_sep.KJ = gamma_03;
%             gamma_sep.KJ2 = gamma_04;
            
        end
        
        function [gamma, res, SH, f_n, mode_shape, k_mesh, gamma_sep] = scaled_sweep(obj_ref, P_scale, n_R_scale, normalize_freq, N_freq)
            %SCALED_SWEEP performs a sweep on the rated power parameter of
            % the Drivetrain object. Returns the scaling factors gamma
            %
            % see also SCALED_VERSION
            %
            
            n_P = numel(P_scale);
            n_fn = numel(obj_ref.modal_analysis);
            
            gamma = zeros(18, n_P);
            res = struct;
            gamma_sep = struct;
            SH = zeros(numel(obj_ref.S_H), n_P);
            f_n = zeros(n_fn, n_P);
            mode_shape = zeros(n_fn, n_fn, n_P);
            k_mesh = zeros(numel(obj_ref.stage), n_P);
            
            gamma_P = P_scale./obj_ref.P_rated;
            
            [~, MS_ref] = obj_ref.modal_analysis;
            
            SH_ref = obj_ref.S_H;
            SS = [SH_ref, SH_ref];
            
            plot_prop1 = {'lineStyle', '-' , 'lineWidth', 2.0, 'color', [346.6667e-3   536.0000e-3   690.6667e-3]};
            plot_prop2 = {'lineStyle', '--', 'lineWidth', 2.0, 'color', [915.2941e-3   281.5686e-3   287.8431e-3]};
            
            figure("units", "centimeters", "position", [5.0 5.0 34.0 12.0]);
            subplot(2, 6, 1:3)
            rectangle(obj_ref);
            xlim([0 6000])
            ylim([-1 1]*1500)
            title(sprintf("Reference: %.1f kW", obj_ref.P_rated));
            
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
%                     bar(SS);
                end
            end
            
            fig_axes = findobj(gcf, "Type", "Axes");
            
            font_setting = {"fontName", "Times", "fontSize", 12.0};
            label_figure(fig_axes, font_setting);
            
            fig_name = @(i)(sprintf("plots\\sweep_scale\\scaled_%d_%d", i, n_P));
            
            file_name = "plots\scale_sweep.gif";
            
            k_P = mean(diff(P_scale));
            
            if(k_P > 0.0)
                gamma_0 = ones(18, 1)*1.0e-3;
            else
                gamma_0 = ones(18, 1);
            end
            
            for idx = 1:n_P
                [obj_sca, gamma_0, res_idx, gamma_idx] = obj_ref.scaled_version(P_scale(idx), n_R_scale, normalize_freq, N_freq, gamma_0);
                
                gamma(:, idx) = gamma_0;
                res(idx).stage = res_idx.stage;
                res(idx).gear  = res_idx.gear;
                res(idx).KJ    = res_idx.KJ;
%                 res(idx).KJ2   = res_idx.KJ2;
                
                gamma_sep(idx).stage = gamma_idx.stage;
                gamma_sep(idx).gear  = gamma_idx.gear;
                gamma_sep(idx).KJ    = gamma_idx.KJ;
%                 gamma_sep(idx).KJ2   = gamma_idx.KJ2;
                
                SH(:, idx) = obj_sca.S_H;
                k_mesh(:, idx) = [obj_sca.stage.k_mesh]';
                [f_n(:, idx), mode_shape(:, :, idx)] = obj_sca.modal_analysis;

                subplot(2, 6, 7:9)
                cla;
                rectangle(obj_sca);
                xlim([0 6000])
                ylim([-1 1]*1500)
                title(sprintf("Scale: %.1f kW = %.2f %% of Ref.", obj_sca.P_rated, gamma_P(idx)*100.0));
                SS(:, 2) = SH(:, idx);
                
                for jdx = 1:3
                    subplot(2, 6, jdx + 3)
                    cla;
                    hold on;
                    plot(1:14, MS_ref(:, jdx)*90.0, plot_prop1{:});
                    plot(1:14, mode_shape(:, jdx, idx)*90.0, plot_prop2{:});

                    if(jdx ~= 3)
                        subplot(2, 6, jdx + 9)
                        cla;
                        hold on;
                        plot(1:14, MS_ref(:, jdx + 3)*90.0, plot_prop1{:})
                        plot(1:14, mode_shape(:, jdx + 3, idx)*90.0, plot_prop2{:});
                    else
                        subplot(2, 6, 12)
                        cla;
                        bar(SS);
                    end
                end

%                 set(fig_axes([4 6]), 'xticklabel', []);
%                 set(fig_axes(2:4), 'yticklabel', []);
                set(fig_axes(2:6), 'xlim', [1 14]);
                set(fig_axes(2:6), 'ylim', [-1 1]*100);
                set(fig_axes(2:6), 'xticklabel', ["C_1", "P_{21}", "S_1", "P_{12}", "P_{32}", "W_3", "G"]);
                
                set(fig_axes, font_setting{:});
                
                savefig(gcf, fig_name(idx));
                print(fig_name(idx), '-dpng');
                saveasGIF(file_name, idx);
                
            end
            
        end
        
        %% Dynamics:
        function H = FRF(obj, freq, varargin)
            
            if(nargin == 2)
                beta = 0.05;
            else
                beta = varargin{1};
            end
            
            switch(obj.dynamic_model)
                case "Thomson_ToV"
                    [M, K_tmp] = obj.Thomson_ToV;
                    
                    K = @(x)(K_tmp);
                    
                    n = length(M);
                    f = ones(n, 1);
                case "Kahraman_1994"
                    [M, K_tmp] = obj.Kahraman_1994;
                    
                    K = @(x)(K_tmp);
                    
                    n = length(M);
                    f = zeros(n, 1);
                    f(1)   = 1.0;
                    f(end) = 1.0;
                case "Lin_Parker_1999"
%                   [M, K, K_b, K_m, K_Omega, G]
                    [M, ~, K_b, K_m, K_Omega, ~] = obj.Lin_Parker_1999;
                    K_tmp = K_b + K_m;
                    
                    K = @(x)(K_b + K_m - K_Omega.*x.^2);
                    
                    n = length(M);
                    f = zeros(n, 1);
                    f(1:3) = 1.0;
                    f(end) = 1.0;
                otherwise
                    error("Option [%s] is NOT available or not implemented yet.", upper(obj.dynamic_model));
                    
%                 case "Eritenel_2011"
%                     [M, K] = obj.E
            end
            
            Omega = 2.0*pi*freq;
            i = sqrt(-1.0);
            
            n_Om = length(Omega);
            H = zeros(n_Om, n);
            
            for idx = 1:n_Om
                H(idx, :) = (-M.*Omega(idx).^2 + i.*Omega(idx).*beta.*K(Omega(idx)) + K(Omega(idx)))\f;
            end
        end
        
        function [f, mode_shape] = nth_resonance(obj, n)
            [f_n, mode_shape] = modal_analysis(obj);
            
            if((n < 1) || (n > numel(f_n)))
                error("n = %d > or ~= 0.");
            end
            
            f = f_n(n);
            mode_shape = mode_shape(:, n);
        end
        
        function [f, mode_shape] = scale_nth_resonance(obj, n, gamma_P, gamma_n, gamma, aspect)
            obj_sca = scale_aspect(obj, gamma_P, gamma_n, gamma, aspect);
            
            [f, mode_shape] = obj_sca.nth_resonance(n);
            
        end
        
        function [f_n, mode_shape] = modal_analysis(obj)
            % MODAL_ANALYSIS calculates the resonances and mode shapes of
            % the Drivetrain via a symmertic eigenvalue problem [1]. The
            % dynamic matrices M and K can be obtained using different
            % modelling approaches, e.g. [2-3].
            %
            % [1] D. Inman, Engineering Vibrations, 4th ed. Boston:
            % Pearson, 2014, pp. 408-413.
            % [2] A. Kahraman, "Natural Modes of Planetary Gear Trains",
            % Journal of Sound and Vibration, vol. 173, no. 1, pp. 125-130,
            % 1994. https://doi.org/10.1006/jsvi.1994.1222.
            % [3] J. Lin and R. Parker, "Analytical Characterization of the
            % Unique Properties of Planetary Gear Free Vibration", Journal
            % of Vibration and Acoustics, vol. 121, no. 3, pp. 316-321,
            % 1999. https://doi.org/10.1115/1.2893982
            %

            switch(obj.dynamic_model)
                case "Thomson_ToV"
                    [M, K] = obj.Thomson_ToV;
                    
                case "Kahraman_1994"
                    [M, K] = obj.Kahraman_1994;
                    
                case "Lin_Parker_1999"
%                   [M, K, K_b, K_m, K_Omega, G]
                    [M, ~, K_b, K_m,       ~, ~] = obj.Lin_Parker_1999;
                    K = K_b + K_m;
                    
                otherwise
                    error("Option [%s] is NOT available or not implemented yet.", upper(obj.dynamic_model));
                    
%                 case "Eritenel_2011"
%                     [M, K] = obj.E
            end
            
            % Symmetric eigenvalue problem [1]:
            L = chol(M, "lower");
            K_tilde = L\K/(L');
            
            % correcting numeric erros and make the problem symmetric:
            K_tilde = (K_tilde + K_tilde')/2.0;

            [mode_shape, D] = eig(K_tilde);
            D = diag(D);   % matrix to vector
            w_n = sqrt(D); % lambda to omega_n

            f_n = w_n./(2.0*pi); % rad/s to Hz
            
            % sorting in ascending order:
            [f_n, idx] = sort(f_n);
            mode_shape = mode_shape(:, idx);
            
            % Normalizing the mode shapes so that the maximum is always +1:
            for idx = 1:length(f_n)
                [ms_max, n] = max(abs(mode_shape(:, idx)));
                mode_shape(:, idx) = mode_shape(:, idx)*sign(mode_shape(n, idx))./ms_max;
            end
            
        end
        
        function [f_n, mode_shape] = scale_modal_analysis(obj, gamma_P, gamma_n, gamma, aspect)
            %SCALED_MODAL_ANALYSIS performs modal analysis on a scaled
            % Drivetrain object, returning only the N first resonances and
            % mode shapes.
            %
            % See also: MODAL_ANALYSIS.
            %
            
            obj_sca = scale_aspect(obj, gamma_P, gamma_n, gamma, aspect);
            
            [f_n, mode_shape] = obj_sca.modal_analysis;
            
        end
        
        function [f_n, mode_shape] = resonances(obj, N, normalize)
            %RESONANCES returns the first N resonances and mode shapes of a
            % Drivetrain object. The resonances can be normalized or not.
            %
            
            [f_n, mode_shape] = obj.modal_analysis;
            
            N_fn = numel(f_n);
            
            if(N <= 0)
                error("N = %d < 0. It should be positive and smaller than %d.", N, N_fn);
            elseif(N > N_fn)
                error("N = %d > %d. It should be positive and smaller than %d.", N, N_fn, N_fn);
            else
                f_n = f_n(1:N);
                mode_shape = mode_shape(:, 1:N);
            end
            
            if(normalize == true)
                f_n = f_n./f_n(1);
                f_n = f_n(2:end);
            end
            
        end
        
        function [f_n, mode_shape] = scale_resonances(obj, N, normalize, gamma_P, gamma_n, gamma, aspect)
            %SCALED_RESONANCES returns the N first resonances and mode
            % shapes of a scaled Drivetrain object. The resonances can be
            % normalized or not.
            %
            
            obj_sca = scale_aspect(obj, gamma_P, gamma_n, gamma, aspect);
            
            [f_n, mode_shape] = obj_sca.resonances(N, normalize);
            
        end
        
        function [M, K] = Thomson_ToV(obj)
            %THOMSON_TOV Returns the inertia and stiffness matrices of the
            % drivetrain according to:
            % W. Thomson and M. Dahleh, Theory of vibration with
            % applications, 5th ed. Prentice-Hall: New Jersey, 1998,
            % pp. 380-381.
            %
            
            J_R = obj.J_Rotor; % [kg-m^2], Rotor inertia
            J_G = obj.J_Gen;   % [kg-m^2], Generator inertia
            
            U = obj.u;
            
            k_LSS = obj.main_shaft.stiffness("torsional");
            k_HSS = obj.stage(end).out_shaft.stiffness("torsional");
            
            k = (k_LSS*k_HSS*U^2)/(k_LSS + k_HSS*U^2);
            
            M = diag([J_R J_G*U^2]);
            K = k*[1.0 -1.0;
                  -1.0  1.0];
        end
        
        function [M, K] = Kahraman_1994(obj)
            %KAHRAMAN_1994 Returns the inertia and stiffness matrices of
            % the drivetrain according to:
            % A. Kahraman, "Natural Modes of Planetary Gear Trains",
            % Journal of Sound and Vibration, vol. 173, no. 1, pp. 125-130,
            % 1994. https://doi.org/10.1006/jsvi.1994.1222.

            J_R = obj.J_Rotor; % [kg-m^2], Rotor inertia according to [3]
            J_G = obj.J_Gen;   % [kg-m^2], Generator inertia according to [4]
            
            M = zeros(14, 14);
            K = zeros(14, 14);
            
            M(1, 1) = J_R;              M(end, end) = J_G;
            
            M_tmp = obj.main_shaft.inertia_matrix("torsional");
            K_tmp = obj.main_shaft.stiffness_matrix("torsional");
            
            M(1:2, 1:2) = M(1:2, 1:2) + M_tmp;        K(1:2, 1:2) = K(1:2, 1:2) + K_tmp;
            
            for idx = 1:2
                [M_tmp, K_tmp] = obj.stage(idx).Kahraman_1994;
                
                jdx = 5*idx - 3;
                kdx = jdx:(jdx + 5);
                M(kdx, kdx) = M(kdx, kdx) + M_tmp;
                K(kdx, kdx) = K(kdx, kdx) + K_tmp;
            end
            
            idx = idx + 1;
            
            [M_tmp, K_tmp] = obj.stage(idx).Kahraman_1994;
            
            jdx = 5*idx - 3;
            kdx = jdx:(jdx + 2);
            M(kdx, kdx) = M(kdx, kdx) + M_tmp;
            K(kdx, kdx) = K(kdx, kdx) + K_tmp;
        end
        
        function [M, K, K_b, K_m, K_Omega, G] = Lin_Parker_1999(obj)
            %LIN_PARKER_1999 Returns the inertia and stiffness matrices of
            % the drivetrain according to:
            % J. Lin and R. Parker, "Analytical Characterization of the
            % Unique Properties of Planetary Gear Free Vibration", Journal
            % of Vibration and Acoustics, vol. 121, no. 3, pp. 316-321,
            % 1999. https://doi.org/10.1115/1.2893982
            %

            m_R = obj.m_Rotor; % [kg], Rotor mass
            m_G = obj.m_Gen;   % [kg], Generator mass according to [?]
            
            J_R = obj.J_Rotor; % [kg-m^2], Rotor inertia according to [3]
            J_G = obj.J_Gen;   % [kg-m^2], Generator inertia according to [4]
            
            M       = zeros(42, 42);
            K_b     = zeros(42, 42);
            K_m     = zeros(42, 42);
            K_Omega = zeros(42, 42);
            G       = zeros(42, 42);
            
            M(1, 1) = m_R;              M(end - 2, end - 2) = m_G;
            M(2, 2) = m_R;              M(end - 1, end - 1) = m_G;
            M(3, 3) = J_R;              M(end    , end    ) = J_G;
            
            M_tmp = obj.main_shaft.inertia_matrix("full");
            K_tmp = obj.main_shaft.stiffness_matrix("full");
            
%             idx = [3 5 6 9 11 12];
            idx = [1 5 6 7 11 12];
            M_tmp(idx, :) = [];                 M_tmp(:, idx) = [];
            K_tmp(idx, :) = [];                 K_tmp(:, idx) = [];
            
%             M_tmp = diag(diag(M_tmp));           K_tmp = diag(diag(K_tmp));
            
            M(  1:6, 1:6) = M(  1:6, 1:6) + M_tmp;
            K_b(1:6, 1:6) = K_b(1:6, 1:6) + K_tmp;
            
            for idx = 1:2
                [M_tmp, ~, K_b_tmp, K_m_tmp, K_Omega_tmp, G_tmp] = obj.stage(idx).Lin_Parker_1999;
                
                jdx = 15*idx - 11;
                kdx = jdx:(jdx + 17);
                
                M(kdx, kdx) = M(kdx, kdx) + M_tmp;
                
                K_b(    kdx, kdx) = K_b(    kdx, kdx) + K_b_tmp;
                K_m(    kdx, kdx) = K_m(    kdx, kdx) + K_m_tmp;
                K_Omega(kdx, kdx) = K_Omega(kdx, kdx) + K_Omega_tmp;
                
                G(kdx, kdx) = G(kdx, kdx) + G_tmp;
                
            end
            
            idx = idx + 1;
            
            [M_tmp, ~, K_b_tmp, K_m_tmp, K_Omega_tmp, G_tmp] = obj.stage(idx).Lin_Parker_1999;
            
            jdx = 15*idx - 11;
            kdx = jdx:(jdx + 8);
            
            M(kdx, kdx) = M(kdx, kdx) + M_tmp;
            
            K_b(    kdx, kdx) = K_b(    kdx, kdx) + K_b_tmp;
            K_m(    kdx, kdx) = K_m(    kdx, kdx) + K_m_tmp;
            K_Omega(kdx, kdx) = K_Omega(kdx, kdx) + K_Omega_tmp;
            
            G(kdx, kdx) = G(kdx, kdx) + G_tmp;
            
            K = @(Om)(K_b + K_m - K_Omega*Om^2);
        end
        
        %% Pitting:
        function [SH_vec, SShaft_vec] = safety_factors(obj)
            %SAFETY_FACTORS returns the safety factors for a Drivetrain
            % object. It calculates the safety factor against pitting
            % for the Gear objects according to ISO 6336-2 and the safety
            % against fatigue and yield according to Shigley's book.
            
            S_Hmin = 1.25;      % [-],  Minimum required safety factor for surface durability according to IEC 61400-4.
            L_h    = 20*365*24; % [h],  Required life
            Q      = 6;         % [-],  ISO accuracy grade
            R_a    = 0.8;       % [um], Maximum arithmetic mean roughness for external gears according to [7], Sec. 7.2.7.2.
            K_A    = 1.25;      % [-],  Application factor
            
            S_ut = Material.S_ut*1.0e-6; % [MPa], Tensile strength
            S_y  = Material.S_y*1.0e-6;  % [MPa], Yield strength
            K_f  = 1.0;                  % [-],   Fatigue stress-concentration factor for bending
            K_fs = 1.0;                  % [-],   Fatigue stress-concentration factor for torsion
            
            T_m  = obj.T_out(1)*obj.stage(1).u;

            SH_vec     = zeros(6, 1);
            SShaft_vec = zeros(4, 1);
            
            SShaft_vec(1) = obj.main_shaft.safety_factors(S_ut, S_y, K_f, K_fs, T_m);
            
            for idx = 1:obj.N_stg
                SShaft_vec(idx + 1) = obj.stage(idx).out_shaft.safety_factors(S_ut, S_y, K_f, K_fs, obj.T_out(idx));
                
                [SH, ~] = obj.stage(idx).Pitting_ISO(obj.P_rated, obj.n_out(idx), S_Hmin, L_h, Q, R_a, K_A);
                
                jdx = 2*idx - 1;
                kdx = jdx:(jdx + 1);
                
                SH_vec(kdx)     = SH;
            end
            
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
            val = zeros(obj.N_stg, 1);
            
            for idx = 1:obj.N_stg
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
        
        function val = get.S_shaft(obj)
            val = obj.S_shaft_val;
        end
    end
    
end
