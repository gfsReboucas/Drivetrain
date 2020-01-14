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
    
    properties(SetAccess = private)
        stage       (1, :) Gear_Set;                                                            % [-],      gearbox stages
        P_rated     (1, :)          {mustBeNumeric, mustBeFinite, mustBePositive} = 5.0e3;      % [kW],     Rated power
        n_rotor     (1, :)          {mustBeNumeric, mustBeFinite, mustBePositive} = 12.1;       % [1/min.], Rated rotor speed
        main_shaft  (1, :) Shaft                                                  = Shaft;      % [-],      Input Shaft
        m_Rotor     (1, :)          {mustBeNumeric, mustBeFinite, mustBePositive} = 110.0e3;    % [kg],     Rotor mass according to [3]
        J_Rotor     (1, :)          {mustBeNumeric, mustBeFinite, mustBePositive} = 57231535.0; % [kg-m^2], Rotor mass moment of inertia according to [6]
        m_Gen       (1, :)          {mustBeNumeric, mustBeFinite, mustBePositive} = 1900.0;     % [kg],     Generator mass according to [4]
        J_Gen       (1, :)          {mustBeNumeric, mustBeFinite, mustBePositive} = 534.116;    % [kg-m^2], Generator mass moment of inertia [4]
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
        T_1;     % [N-m],    Output torque for each stage
        n_1;     % [1/min.], Output speed  for each stage
        u;       % [-],      Cumulative gear ratio
        S_H;     % [-],      Safety factor for surface durability (against pitting)
        S_shaft; % [-],      Safey factor for the shafts
    end
    
    methods
        function obj = Drivetrain(stage, P_r, n_r, inp_shaft, m_R, J_R, m_G, J_G)
            if(nargin == 0)
                for idx = 1:3
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
            
            [obj.S_H_val, obj.S_shaft_val] = obj.safety_factors;
            
        end
        
        function tab = disp(obj)
            % description, symbol, unit, value
            
            tab_str = {"Rated power",                                  "P",       "kW",     "-+-",                obj.P_rated,          "-+-";            % 1
                       "Speed (Sun/Pinion)",                           "n_1",     "1/min.", obj.n_1(1),           obj.n_1(2),           obj.n_1(3);
                       "Torque (Sun/Pinion)",                          "T_1",     "N-m",    obj.T_1(1),           obj.T_1(2),           obj.T_1(3);
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
                fprintf("Drivetrain:\n");
                disp(tab);
                fprintf("Main shaft:\n");
                disp(obj.main_shaft);
                clear tab;
            end

        end
        
        function [tab, tab_str] = comp_stage(ref, sca, idx)
            tab_str = {"Rated power",                                  "P",       "kW",     ref.P_rated,                sca.P_rated,                sca.P_rated               /ref.P_rated;                % 1
                       "Speed (Sun/Pinion)",                           "n_1",     "1/min.", ref.n_1(idx),               sca.n_1(idx),               sca.n_1(idx)              /ref.n_1(idx);               % 2
                       "Torque (Sun/Pinion)",                          "T_1",     "N-m",    ref.T_1(idx),               sca.T_1(idx),               sca.T_1(idx)              /ref.T_1(idx);               % 3
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
            for idx = 1:3
                subplot(1, 3, idx)
                obj.stage(idx).plot;
            end
        end
        
        function rectangle(obj)
            addpath("\\home.ansatt.ntnu.no\geraldod\Documents\MATLAB\Plot\linspecer");
            color = linspecer(6, "qualitative");
            
            hold on;
            C_s = [obj.main_shaft.L/2.0 0.0]';
            rectangle(obj.main_shaft, C_s, color(6, :), "edgeColor", "k", "lineStyle", "-" , "faceColor", color(6, :));
            C = [obj.main_shaft.L;
                 obj.stage(1).carrier.b + obj.stage(1).b + obj.stage(1).out_shaft.L;
                 obj.stage(2).carrier.b + obj.stage(2).b + obj.stage(2).out_shaft.L];

            for idx = 1:3
                rectangle(obj.stage(idx), [sum(C(1:idx)) 0.0]');
            end
            hold off;
        end
        
        function plot_comp(DT1, DT2)
            hold on;
            axis equal;
            box on;
            for idx = 1:3
                subplot(2, 3, idx)
                DT1.stage(idx).plot;
                title(sprintf("Stage %d", idx));
                subplot(2, 3, idx + 3)
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
            
            if(numel(gamma) < 18)
                error("gamma must have 18 elements.");
            end
            
            gamma_s1_mn = gamma( 1, :); % Stage 01, normal module
            gamma_s1_bb = gamma( 2, :); % Stage 01, face width
            gamma_s1_LL = gamma( 3, :); % Stage 01, output shaft length
%             gamma_s1_dd = gamma( 4, :); % Stage 01, output shaft diameter
            gamma_s2_mn = gamma( 5, :); % Stage 02, normal module
            gamma_s2_bb = gamma( 6, :); % Stage 02, face width
            gamma_s2_LL = gamma( 7, :); % Stage 02, output shaft length
%             gamma_s2_dd = gamma( 8, :); % Stage 02, output shaft diameter
            gamma_s3_mn = gamma( 9, :); % Stage 03, normal module
            gamma_s3_bb = gamma(10, :); % Stage 03, face width
            gamma_s3_LL = gamma(11, :); % Stage 02, output shaft length
%             gamma_s3_dd = gamma(12, :); % Stage 03, output shaft diameter
            gamma_R_m   = gamma(13, :); % Rotor, mass
            gamma_R_J   = gamma(14, :); % Rotor, mass moment of inertia
            gamma_G_m   = gamma(15, :); % Generator, mass
            gamma_G_J   = gamma(16, :); % Generator, mass moment of inertia
%             gamma_S_L     = gamma(17, :); % Main shaft, length
%             gamma_S_d     = gamma(18, :); % Main shaft, diameter
            
            gamma_Torque = gamma_P/gamma_n; % Applied torque
            
            % Shaft diameters:
            gamma_shaft_d = nthroot(gamma_Torque, 3.0);
            
            gamma_s1_dd = gamma_shaft_d;
            gamma_s2_dd = gamma_shaft_d;
            gamma_s3_dd = gamma_shaft_d;
            
            gamma_S_d = gamma_s3_dd; % Main shaft, diameter
            gamma_S_L = gamma_s3_LL; % Main shaft, length
            
            main_shaft_sca = Shaft(obj_ref.main_shaft.d*gamma_S_d, ...
                                   obj_ref.main_shaft.L*gamma_S_L);
            
            stage_sca = [Gear_Set, Gear_Set, Gear_Set];
            
            gamma_mn = [gamma_s1_mn gamma_s2_mn gamma_s3_mn]';
            gamma_bb = [gamma_s1_bb gamma_s2_bb gamma_s3_bb]';
            gamma_LL = [gamma_s1_LL gamma_s2_LL gamma_s3_LL]';
            gamma_dd = [gamma_s1_dd gamma_s2_dd gamma_s3_dd]';
            
            for idx = 1:3
                % Scaled stage idx:
                gamma_stage = [gamma_mn(idx) gamma_bb(idx)]';
                % Scaled output shaft idx:
                gamma_shaft = [gamma_dd(idx) gamma_LL(idx)]';
                stage_sca(idx) = obj_ref.stage(idx).scaled_Gear_Set(gamma_stage, gamma_shaft);
            end
            
            obj_sca = Drivetrain(stage_sca, ...
                                 obj_ref.P_rated *gamma_P, ...
                                 obj_ref.n_rotor *gamma_n, ...
                                 main_shaft_sca, ...
                                 obj_ref.m_Rotor*gamma_R_m, ...
                                 obj_ref.J_Rotor*gamma_R_J, ...
                                 obj_ref.m_Gen  *gamma_G_m, ...
                                 obj_ref.J_Gen  *gamma_G_J);
        end
        
        function obj_sca = scale_by_stage(obj_ref, gamma_P, gamma_n, gamma)
            %SCALE_BY_STAGE returns a Drivetrain object scaled by the
            % factors gamma_P and gamma_n for its rated power and rotor
            % speed and is scaled by stage (i.e. one scale factor for both
            % the normal module and face width per stage).
            %
            
            if(numel(gamma) < 3)
                error("gamma must have 3 elements.");
            end
            
            gamma_s1 = gamma(1, :);
            gamma_s2 = gamma(2, :);
            gamma_s3 = gamma(3, :);
            
            gamma_full = [gamma_s1, gamma_s1, 1.0, 1.0, ...
                          gamma_s2, gamma_s2, 1.0, 1.0, ...
                          gamma_s3, gamma_s3, 1.0, 1.0, ...
                          ones(1, 6)]';

            obj_sca = scale_Drivetrain(obj_ref, gamma_P, gamma_n, gamma_full);
            
        end
        
        function obj_sca = scale_gear(obj_ref, gamma_P, gamma_n, gamma)
            %SCALE_GEAR returns a Drivetrain object scaled by the factors
            % gamma_P and gamma_n for its rated power and rotor speed and
            % scales only the gear dimensions (normal module and face
            % width).
            %
            
            if(numel(gamma) < 6)
                error("gamma must have 6 elements.");
            end
            
            gamma_s1_mn = gamma(1, :);           gamma_s1_bb = gamma(2, :);
            gamma_s2_mn = gamma(3, :);           gamma_s2_bb = gamma(4, :);
            gamma_s3_mn = gamma(5, :);           gamma_s3_bb = gamma(6, :);
            
            gamma_full = [gamma_s1_mn, gamma_s1_bb, 1.0, 1.0, ...
                          gamma_s2_mn, gamma_s2_bb, 1.0, 1.0, ...
                          gamma_s3_mn, gamma_s3_bb, 1.0, 1.0, ...
                          ones(1, 6)]';

            obj_sca = scale_Drivetrain(obj_ref, gamma_P, gamma_n, gamma_full);
            
        end
        
        function obj_sca = scale_stiffness(obj_ref, gamma_P, gamma_n, gamma)
            %SCALE_STIFFNESS returns a Drivetrain object scaled by the
            % factors gamma_P and gamma_n for its rated power and rotor
            % speed and scales only stiffness-related parameters. The mesh
            % stiffness is assumed to be proportinal to the gear's face
            % width.
            %
            
            if(numel(gamma) < 6)
                error("gamma must have 6 elements.");
            end
            
            gamma_b1 = gamma(1, :);            gamma_L1 = gamma(2, :);
            gamma_b2 = gamma(3, :);            gamma_L2 = gamma(4, :);
            gamma_b3 = gamma(5, :);            gamma_L3 = gamma(6, :);
            
            gamma_full = [1.0, gamma_b1, gamma_L1, 1.0, ...
                          1.0, gamma_b2, gamma_L2, 1.0, ...
                          1.0, gamma_b3, gamma_L3, 1.0, ...
                          ones(1, 6)]';

            obj_sca = scale_Drivetrain(obj_ref, gamma_P, gamma_n, gamma_full);
            
        end
                
        function obj_sca = scale_shaft_stiffness(obj_ref, gamma_P, gamma_n, gamma)
            %SCALE_SHAFT_STIFFNESS returns a Drivetrain object scaled by 
            % the factors gamma_P and gamma_n for its rated power and rotor
            % speed and scales only parameters related to the shaft's 
            % stiffness.
            %
            
            if(numel(gamma) < 3)
                error("gamma must have 3 elements.");
            end
            
            gamma_s1 = gamma(1, :);
            gamma_s2 = gamma(2, :);
            gamma_s3 = gamma(3, :);
            
            gamma_full = [1.0, 1.0, gamma_s1, 1.0, ...
                          1.0, 1.0, gamma_s2, 1.0, ...
                          1.0, 1.0, gamma_s3, 1.0, ...
                          ones(1, 6)]';

            obj_sca = scale_Drivetrain(obj_ref, gamma_P, gamma_n, gamma_full);
            
        end
        
        function obj_sca = scale_shaft_stiffness_mass_mom_inertia(obj_ref, gamma_P, gamma_n, gamma)
            %SCALE_SHAFT_STIFFNESS_MASS_MOM_INERTIA returns a Drivetrain
            % object scaled by the factors gamma_P and gamma_n for its 
            % rated power and rotor speed and scales only parameters 
            % related to the shaft's stiffness and mass momen of inertia of
            % rotor and generator.
            %
            
            if(numel(gamma) < 5)
                error("gamma must have 5 elements.");
            end
            
            gamma_s1 = gamma(1, :);
            gamma_s2 = gamma(2, :);
            gamma_s3 = gamma(3, :);
            gamma_JR = gamma(4, :);
            gamma_JG = gamma(5, :);
            
            gamma_full = [1.0, 1.0     , gamma_s1, 1.0, ...
                          1.0, 1.0     , gamma_s2, 1.0, ...
                          1.0, 1.0     , gamma_s3, 1.0, ...
                          1.0, gamma_JR, 1.0     , gamma_JG, ...
                          1.0, 1.0]';

            obj_sca = scale_Drivetrain(obj_ref, gamma_P, gamma_n, gamma_full);
            
        end
        
        function obj_sca = scale_aspect(obj_ref, gamma_P, gamma_n, gamma, aspect)
            %SCALE_ASPECT returns a Drivetrain object scaled by the factors
            % gamma_P and gamma_n for its rated power and rotor speed. The
            % argument gamma can be used to scale the Drivetrain
            % considering various aspects.
            %
            
            switch(aspect)
                case "Drivetrain"
                    obj_sca = obj_ref.scale_Drivetrain(gamma_P, gamma_n, gamma);
                case "by_stage"
                    obj_sca = obj_ref.scale_by_stage(gamma_P, gamma_n, gamma);
                case "gear"
                    obj_sca = obj_ref.scale_gear(gamma_P, gamma_n, gamma);
                case "stiffness_mass_mom_inertia"
                    obj_sca = obj_ref.scale_shaft_stiffness_mass_mom_inertia(gamma_P, gamma_n, gamma);
                case "shaft_stiffness"
                    obj_sca = obj_ref.scale_shaft_stiffness(gamma_P, gamma_n, gamma);
                case "stiffness"
                    obj_sca = obj_ref.scale_stiffness(gamma_P, gamma_n, gamma);
                otherwise
                    error("Option [%s] is NOT valid.", aspect);
            end
        end
           
        function [obj_sca, gamma, res] = scaled_version(obj_ref, P_scale, n_R_scale, normalize_freq, N_freq)
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
            
            % Input scaling factors:
            gamma_P = P_scale/obj_ref.P_rated;
            gamma_n = n_R_scale/obj_ref.n_rotor;
            
            %% 1. Stage scaling:
            aspect = "by_stage";
            
            S_H_ref = obj_ref.S_H;
            
            fun_SH  = @(x)(S_H_ref - obj_ref.scaled_safety_factors(gamma_P, gamma_n, x, aspect));
            fun_eq = @(x)(norm(fun_SH(x))^2);
            
            gamma_min = ones(3, 1)*1.0e-6;
            gamma_Max = ones(3, 1);
            gamma_0 = (0*gamma_min + 2*gamma_Max)/2.0;
            
            fun_ineq = @(x)(1.25 - obj_ref.scaled_safety_factors(gamma_P, gamma_n, x, aspect));

            constraint_fun = @(x)deal(fun_ineq(x), fun_SH(x));
            
            opt_solver = optimoptions("fmincon", "display", "notify");
            
            id_1 = "prog:input";
            id_2 = "MATLAB:nearlySingularMatrix";
            warning("off", id_1);
            warning("off", id_2);
            [gamma_stage, res_stage] = fmincon(fun_eq, gamma_0, [], [], [], [], gamma_min, gamma_Max, constraint_fun, opt_solver);
            
            %% 2. Gears:
            aspect = "gear";
            
            A = zeros(6, 3);
            A(1, 1) = 1.0;      A(2, 1) = 1.0;
            A(3, 2) = 1.0;      A(4, 2) = 1.0;
            A(5, 3) = 1.0;      A(6, 3) = 1.0;
            gamma_0 = A*gamma_stage;
            
            gamma_min = ones(6, 1)*1.0e-6;
            gamma_Max = ones(6, 1);
            
            fun_SH  = @(x)(S_H_ref - obj_ref.scaled_safety_factors(gamma_P, gamma_n, x, aspect));
            fun_eq = @(x)(norm(fun_SH(x))^2);
            
            fun_ineq = @(x)(1.25 - obj_ref.scaled_safety_factors(gamma_P, gamma_n, x, aspect));

            constraint_fun = @(x)deal(fun_ineq(x), fun_SH(x));
            
            [gamma_gear, res_gear] = fmincon(fun_eq, gamma_0, [], [], [], [], gamma_min, gamma_Max, constraint_fun, opt_solver);
            warning("on", id_1);
            warning("on", id_2);
            
            %% 3. Shaft and Mass moment of inertia:
            f_n_ref = obj_ref.resonances(N_freq, normalize_freq);
            aspect = "stiffness_mass_mom_inertia";
            
            fun_fn = @(x)(1.0 - obj_ref.scale);
            
        end
        
        %% Dynamics:
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

            f_n = w_n'./(2.0*pi); % rad/s to Hz
            
            % sorting in ascending order:
            [f_n, idx] = sort(f_n);
            mode_shape = mode_shape(:, idx);
            
            % Normalizing the mode shapes so that the maximum is always +1:
            for idx = 1:length(f_n)
                [ms_max, n] = max(abs(mode_shape(:, idx)));
                mode_shape(:, idx) = mode_shape(:, idx)*sign(mode_shape(n, idx))./ms_max;
            end
            
        end
        
        function [f_n, mode_shape] = scaled_modal_analysis(obj, gamma_P, gamma_n, gamma, aspect)
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
%                 f_n = f_n(2:end);
            end
            
        end
        
        function [f_n, mode_shape] = scaled_resonances(obj, N, normalize, gamma_P, gamma_n, gamma)
            %SCALED_RESONANCES returns the N first resonances and mode
            % shapes of a scaled Drivetrain object. The resonances can be
            % normalized or not.
            %
            
            [f_n, mode_shape] = scaled_modal_analysis(obj, gamma_P, gamma_n, gamma);
            
            N_fn = numel(f_n);
            
            if(N <= 0)
                error("N = %d < 0. It should be positive and smaller than %d.", N, N_fn);
            elseif(N > N_fn)
                error("N = %d > %d. It should be positive and smaller than %d.", N, N_fn, N_fn);
            else
                f_n = f_n(1:N);
            end
            
            if(normalize == true)
                f_n = f_n./f_n(1);
%                 f_n = f_n(2:end);
            end

        end
        
        function f_n = natural_freq(obj, calc_method, opt_freq, N, gamma_d, gamma)
            
            for idx = 1:3
                sca_shaft = obj.stage(idx).out_shaft.scaled_shaft([gamma_d, gamma(idx)]);
                obj.stage(idx).out_shaft = sca_shaft;
            end
            
            LSS_sca = obj.main_shaft.scaled_shaft([gamma_d, gamma(3)]);
            obj.main_shaft = LSS_sca;
            
            if(strcmp(calc_method, "Kahraman_1994"))
                gm_L = gamma(4:5);
            elseif(strcmp(calc_method, "Lin_Parker_1999"))
                gm_L = gamma(4:7);
            end
            
            f_n = obj.modal_analysis(calc_method, opt_freq, N, gm_L);
        end
        
        function [f_n, eig_vec] = modal_analysis_OLD(obj, calc_method, opt_freq, n_f, gamma)
            %MODAL_ANALYSIS calculates the first n_f resonances and mode
            % shapes of the drivetrain via an symmetric eigenvalue problem
            % [1]. The dynamic matrices M and K are obtained using
            % different modelling approaches, e.g. [2-3].
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
            
            switch(calc_method)
                case "Thomson_ToV"
                    if(~exist("gamma", "var"))
                        gamma = ones(1, 2);
                    end
                    
                    [M, K] = obj.Thomson_ToV(gamma, "");
                case "Kahraman_1994"
                    if(~exist("gamma", "var"))
                        gamma = ones(1, 2);
                    end
                    
                    [M, K] = obj.Kahraman_1994(gamma);
                case "Lin_Parker_1999"
                    if(~exist("gamma", "var"))
                        gamma = ones(1, 4);
                    end
                    
%                   [M, K, K_b, K_m, K_Omega, G]
                    [M, ~, K_b, K_m,       ~, ~] = obj.Lin_Parker_1999(gamma);
                    K = K_b + K_m;
                case "Eritenel_2011"
                    if(~exist("gamma", "var"))
                        gamma = ones(1, 4);
                    end
                    
                    [M, K] = obj.Eritenel_2011(gamma);
                otherwise
                    error("Option [%s] is NOT available.", upper(calc_method));
            end
            
            % Symmetric eigenvalue problem according to Inman:
            L = chol(M, "lower");
            K_tilde = L\K/(L');
            
            K_tilde = (K_tilde + K_tilde')/2.0;

            [eig_vec, D] = eig(K_tilde);
            D = diag(D);   % matrix to vector
            w_n = sqrt(D); % lambda to omega_n

            f_n = w_n'./(2.0*pi); % rad/s to Hz
            
            [f_n, idx] = sort(f_n); % ascending sorting
            eig_vec = eig_vec(:, idx);

            
            if(any(imag(f_n) ~= 0))
                idx_im =  find(imag(f_n) ~= 0);
                idx = 1:length(f_n);
                idx(idx_im) = [];
                idx_2 = [idx, idx_im];
                f_n = f_n(idx_2);
                eig_vec = eig_vec(:, idx_2);
            end
            
            % Normalizing eigenvectors:
            for idx = 1:length(f_n)
                [ev_max, n] = max(abs(eig_vec(:, idx)));
                eig_vec(:, idx) = eig_vec(:, idx)*sign(eig_vec(n, idx))./ev_max;
            end
            
            if(n_f < length(f_n))
                f_n = f_n(1:n_f);
%                 eig_vec = eig_vec(:,1:n_f);
            end
            
            if(strcmp(opt_freq, "freq_ratio"))
                f_n = f_n./f_n(1);
                f_n = f_n(2:end);
            end
            
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
            k_HSS = obj.stage(3).out_shaft.stiffness("torsional");
            
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
            
            T_m  = obj.T_1(1)*obj.stage(1).u;

            SH_vec     = zeros(6, 1);
            SShaft_vec = zeros(4, 1);
            
            SShaft_vec(1) = obj.main_shaft.safety_factors(S_ut, S_y, K_f, K_fs, T_m);
            
            for idx = 1:3
                SShaft_vec(idx + 1) = obj.stage(idx).out_shaft.safety_factors(S_ut, S_y, K_f, K_fs, obj.T_1(idx));
                
                [SH, ~] = obj.stage(idx).Pitting_ISO(obj.P_rated, obj.n_1(idx), S_Hmin, L_h, Q, R_a, K_A);
                
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
    
    %% Get methods:
    methods
        function val = get.n_1(obj)
            val = zeros(1, 3);
            
            for idx = 1:3
                val(idx) = obj.n_rotor*prod([obj.stage(1:idx).u]);
            end
        end
        
        function val = get.T_1(obj)
            val = (obj.P_rated*1.0e3)./(obj.n_1*pi/30.0);
        end
        
        function val = get.u(obj)
            val = prod([obj.stage(1:end).u]);
        end
        
        function val = get.S_H(obj)
            val = obj.S_H_val;
        end
        
        function val = get.S_shaft(obj)
            val = obj.S_shaft_val;
        end
    end
    
end
