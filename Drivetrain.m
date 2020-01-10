classdef Drivetrain
    % This class implements SOME procedures for the dynamic analysis and
    % scaling of drivetrains. The safety factor for surface durability
    % (pitting) is calculated according to ISO 6336 [1, 2]. The NREL 5MW
    % reference gearbox proposed by Nejad et. al. [3] is used as the
    % default drivetrain model, but other models should be implemented in
    % the future.
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
        stage      (1, 3) Gear_Set;                                                                                                     % [-],      gearbox stages
        P_rated    (1, :)          {mustBeNumeric, mustBeFinite, mustBePositive}                                          = 5.0e3;      % [kW],     Rated power
        n_rotor    (1, :)          {mustBeNumeric, mustBeFinite, mustBePositive}                                          = 12.1;       % [1/min.], Rated rotor speed
        input_shaft(1, :) Shaft                                                                                           = Shaft;      % [-],      Input Shaft
        S_shaft    (1, 4)          {mustBeNumeric, mustBeFinite, mustBePositive, mustBeGreaterThanOrEqual(S_shaft, 1.0)}  = 1.0;        % [-],      Safey factor for the shafts
        m_Rotor    (1, :)          {mustBeNumeric, mustBeFinite, mustBePositive}                                          = 110.0e3;    % [kg],     Rotor mass according to [3]
        J_Rotor    (1, :)          {mustBeNumeric, mustBeFinite, mustBePositive}                                          = 57231535.0; % [kg-m^2], Rotor mass moment of inertia according to [6]
        m_Gen      (1, :)          {mustBeNumeric, mustBeFinite, mustBePositive}                                          = 1900.0;     % [kg],     Generator mass according to [4]
        J_Gen      (1, :)          {mustBeNumeric, mustBeFinite, mustBePositive}                                          = 534.116;    % [kg-m^2], Generator mass moment of inertia [4]
%         S_H        (1, 6)          {mustBeNumeric, mustBeFinite, mustBePositive}                                          = 1.25;       % [-],      Safety factor for surface durability (against pitting)
%         sigma_H    (1, 6)          {mustBeNumeric, mustBeFinite, mustBePositive}                                          = 1500.0;     % [N/mm^2], Contact stress
        S_H        (1, :)          {mustBeNumeric, mustBeFinite, mustBePositive, mustBeGreaterThanOrEqual(S_H, 1.0)}      = 1.25;   % [-],      Safety factor for surface durability (against pitting)
        sigma_H    (1, :)          {mustBeNumeric, mustBeFinite, mustBePositive, mustBeLessThanOrEqual(sigma_H, 1500.0)}  = 1500.0; % [N/mm^2], Contact stress
    end
    
    properties(Dependent)
        T_1; % [N-m],    Output torque for each stage
        n_1; % [1/min.], Output speed  for each stage
        u;   % [-],      Cumulative gear ratio
    end
    
    methods
        function obj = Drivetrain(stage, P_r, n_r, inp_shaft)
            if(nargin == 0)
                for idx = 1:3
                    stage(idx) = Gear_Set.NREL_5MW(idx);
                end
                
                P_r = 5.0e3; % [kW], Rated power
                n_r = 12.1; % [1/min.], Input speed
                inp_shaft = Shaft;
            end
            
            obj.stage = stage;
            obj.P_rated = P_r;
            obj.n_rotor = n_r;
            obj.input_shaft = inp_shaft;
            
            S_Hmin = 1.25;      % [-], Minimum required safety factor for surface durability according to 
            L_h    = 20*365*24; % [h],  Required life
            Q      = 6;         % [-],  ISO accuracy grade
            R_a    = 0.8;       % [um], Maximum arithmetic mean roughness for external gears according to [7], Sec. 7.2.7.2.
            K_A    = 1.25;      % [-], Application factor
            
            S_ut = Material.S_ut*1.0e-6; % [MPa], Tensile strength
            S_y  = Material.S_y*1.0e-6;  % [MPa], Yield strength
            K_f  = 1.0;   % [-],   Fatigue stress-concentration factor for bending
            K_fs = 1.0;   % [-],   Fatigue stress-concentration factor for torsion
            T_m  = obj.T_1(1)*obj.stage(1).u;

            obj.S_shaft = obj.input_shaft.safety_factors(S_ut, S_y, K_f, K_fs, T_m);
            
            for idx = 1:3
                obj.S_shaft(idx + 1) = obj.stage(idx).shaft.safety_factors(S_ut, S_y, K_f, K_fs, obj.T_1(idx));
                
                [SH, sigmaH] = obj.stage(idx).Pitting_ISO(obj.P_rated, obj.n_1(idx), S_Hmin, L_h, Q, R_a, K_A);
                
                jdx = 2*idx - 1;
                kdx = jdx:(jdx + 1);
                obj.S_H(kdx) = SH;
                obj.sigma_H(kdx) = sigmaH;
            end
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
                fprintf("Input shaft:\n");
                disp(obj.input_shaft);
                clear tab;
            end

        end
        
        function [tab, tab_str] = comp_stage(ref, sca, idx)
            tab_str = {"Rated power",                                  "P",       "kW",     ref.P_rated,            sca.P_rated,            sca.P_rated           /ref.P_rated;            % 1
                       "Speed (Sun/Pinion)",                           "n_1",     "1/min.", ref.n_1(idx),           sca.n_1(idx),           sca.n_1(idx)          /ref.n_1(idx);           % 2
                       "Torque (Sun/Pinion)",                          "T_1",     "N-m",    ref.T_1(idx),           sca.T_1(idx),           sca.T_1(idx)          /ref.T_1(idx);           % 3
                       "Safety factor against pitting (Sun/Pinion)",   "S_H1",    "-",      ref.S_H(2*idx - 1),     sca.S_H(2*idx - 1),     sca.S_H(2*idx - 1)    /ref.S_H(2*idx - 1);     % 4
                       "Safety factor against pitting (Planet/Wheel)", "S_H2",    "-",      ref.S_H(2*idx),         sca.S_H(2*idx),         sca.S_H(2*idx)        /ref.S_H(2*idx);         % 5
                       "Gear ratio",                                   "u",       "-",      ref.stage(idx).u,       sca.stage(idx).u,       sca.stage(idx).u      /ref.stage(idx).u;       % 6
                       "Normal module",                                "m_n",     "mm",     ref.stage(idx).m_n,     sca.stage(idx).m_n,     sca.stage(idx).m_n    /ref.stage(idx).m_n;     % 7
                       "Face width",                                   "b",       "mm",     ref.stage(idx).b,       sca.stage(idx).b,       sca.stage(idx).b      /ref.stage(idx).b;       % 8
                       "Center distance",                              "a_w",     "mm",     ref.stage(idx).a_w,     sca.stage(idx).a_w,     sca.stage(idx).a_w    /ref.stage(idx).a_w;     % 9
                       "Reference diameter (Sun/Pinion)",              "d_1",     "mm",     ref.stage(idx).d(1),    sca.stage(idx).d(1),    sca.stage(idx).d(1)   /ref.stage(idx).d(1);    % 10
                       "Mass (Sun/Pinion)",                            "m_1",     "kg",     ref.stage(idx).mass(1), sca.stage(idx).mass(1), sca.stage(idx).mass(1)/ref.stage(idx).mass(1); % 11
                       "Mass moment of inertia (Sun/Pinion)",          "J_xx1",   "kg-m^2", ref.stage(idx).J_x(1),  sca.stage(idx).J_x(1),  sca.stage(idx).J_x(1) /ref.stage(idx).J_x(1);  % 12
                       "Shaft / Diameter",                             "d",       "mm",     ref.stage(idx).shaft.d, sca.stage(idx).shaft.d, sca.stage(idx).shaft.d/ref.stage(idx).shaft.d; % 13
                       "Shaft / Length",                               "L",       "mm",     ref.stage(idx).shaft.L, sca.stage(idx).shaft.L, sca.stage(idx).shaft.L/ref.stage(idx).shaft.L; % 14
                        };
%                        "Reference diameter (Planet/Wheel)",            "d_2",     "mm",     ref.stage(idx).d(2),    sca.stage(idx).d(2),    sca.stage(idx).d(2)   /ref.stage(idx).d(2);    % 11
%                        "Reference diameter (Ring)",                    "d_3",     "mm",     "-*-",                  "-*-",                  "-*-";                                         % 12
%                        "Mass (Planet/Wheel)",                          "m_2",     "kg",     ref.stage(idx).mass(2), sca.stage(idx).mass(2), sca.stage(idx).mass(2)/ref.stage(idx).mass(2); % 14
%                        "Mass (Ring)",                                  "m_3",     "kg",     "-*-",                  "-*-",                  "-*-";                                         % 15
%                        "Mass moment of inertia (Planet/Wheel)",        "J_xx2",   "kg-m^2", ref.stage(idx).J_x(2),  sca.stage(idx).J_x(2),  sca.stage(idx).J_x(2) /ref.stage(idx).J_x(2);  % 17
%                        "Mass moment of inertia (Ring)",                "J_xx3",   "kg-m^2", "-*-",                  "-*-",                  "-*-";                                         % 18

                    
%             if(idx < 3)
%                 tab_str(12, :) = {"Reference diameter (Ring)",                    "d_3",     "mm",     ref.stage(idx).d(3),    sca.stage(idx).d(3),    sca.stage(idx).d(3)   /ref.stage(idx).d(3)};
%                 tab_str(15, :) = {"Mass (Ring)",                                  "m_3",     "kg",     ref.stage(idx).mass(3), sca.stage(idx).mass(3), sca.stage(idx).mass(3)/ref.stage(idx).mass(3)};
%                 tab_str(18, :) = {"Mass moment of inertia (Ring)",                "J_xx3",   "kg-m^2", ref.stage(idx).J_x(3),  sca.stage(idx).J_x(3),  sca.stage(idx).J_x(3) /ref.stage(idx).J_x(3)};
%             end
            
            Reference = tab_str(:,4);
            Scale     = tab_str(:,5);
            Ratio     = tab_str(:,6);
            
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
            C_s = [obj.input_shaft.L/2.0 0.0]';
            rectangle(obj.input_shaft, C_s, color(6, :), "edgeColor", "k", "lineStyle", "-" , "faceColor", color(6, :));
            C = [obj.input_shaft.L;
                 obj.stage(1).carrier.b + obj.stage(1).b + obj.stage(1).shaft.L;
                 obj.stage(2).carrier.b + obj.stage(2).b + obj.stage(2).shaft.L];

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
        function [obj_scale, gamma, res] = scale(obj_ref, P_scale)
            gamma_P = P_scale/obj_ref.P_rated;
            
        end
        
        %% Dynamics:
        function f_n = natural_freq(obj, calc_method, opt_freq, N, gamma_d, gamma)
            
            for idx = 1:3
                sca_shaft = obj.stage(idx).shaft.scaled_shaft([gamma_d, gamma(idx)]);
                obj.stage(idx).shaft = sca_shaft;
            end
            
            LSS_sca = obj.input_shaft.scaled_shaft([gamma_d, gamma(3)]);
            obj.input_shaft = LSS_sca;
            
            if(strcmp(calc_method, "Kahraman_1994"))
                gm_L = gamma(4:5);
            elseif(strcmp(calc_method, "Lin_Parker_1999"))
                gm_L = gamma(4:7);
            end
            
            f_n = obj.modal_analysis(calc_method, opt_freq, N, gm_L);
        end
        
        function [f_n, eig_vec] = modal_analysis(obj, calc_method, opt_freq, n_f, gamma)
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
        
        function [M, K] = Thomson_ToV(obj, gamma, brake)
            %THOMSON_TOV Returns the inertia and stiffness matrices of the
            % drivetrain according to:
            % W. Thomson and M. Dahleh, Theory of vibration with
            % applications, 5th ed. Prentice-Hall: New Jersey, 1998,
            % pp. 380-381.
            %

            if(~exist("gamma", "var"))
                gamma_J_R = 1.0;
                gamma_J_G = 1.0;
            else
                gamma_J_R = gamma(1);
                gamma_J_G = gamma(2);
            end
            
            J_R = obj.J_Rotor; % [kg-m^2], Rotor inertia according to [3]
            J_G = obj.J_Gen;   % [kg-m^2], Generator inertia according to [4]
            
            J_R = J_R*gamma_J_R;
            J_G = J_G*gamma_J_G;
            
            U = obj.u;
            
            k_LSS = obj.input_shaft.stiffness("torsional");
            k_HSS = obj.stage(3).shaft.stiffness("torsional");
            
            k = (k_LSS*k_HSS*U^2)/(k_LSS + k_HSS*U^2);
            
            M = diag([J_R J_G*U^2]);
            K = k*[1.0 -1.0;
                  -1.0  1.0];
              
            if(strcmp(brake, "brake"))
                M = M(1, 1);
                K = K(1, 1);
            end
            
        end
        
        function [M, K] = Kahraman_1994(obj, gamma)
            %KAHRAMAN_1994 Returns the inertia and stiffness matrices of
            % the drivetrain according to:
            % A. Kahraman, "Natural Modes of Planetary Gear Trains",
            % Journal of Sound and Vibration, vol. 173, no. 1, pp. 125-130,
            % 1994. https://doi.org/10.1006/jsvi.1994.1222.

            if(~exist("gamma", "var"))
                gamma_J_R = 1.0;
                gamma_J_G = 1.0;
            else
                gamma_J_R = gamma(1);
                gamma_J_G = gamma(2);
            end
            
            J_R = obj.J_Rotor; % [kg-m^2], Rotor inertia according to [3]
            J_G = obj.J_Gen;   % [kg-m^2], Generator inertia according to [4]
            
            J_R = J_R*gamma_J_R;
            J_G = J_G*gamma_J_G;
            
            M = zeros(14, 14);
            K = zeros(14, 14);
            
            M(1, 1) = J_R;              M(end, end) = J_G;
            
            M_tmp = obj.input_shaft.inertia_matrix("torsional");
            K_tmp = obj.input_shaft.stiffness_matrix("torsional");
            
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
        
        function [M, K, K_b, K_m, K_Omega, G] = Lin_Parker_1999(obj, gamma)
            %LIN_PARKER_1999 Returns the inertia and stiffness matrices of
            % the drivetrain according to:
            % J. Lin and R. Parker, "Analytical Characterization of the
            % Unique Properties of Planetary Gear Free Vibration", Journal
            % of Vibration and Acoustics, vol. 121, no. 3, pp. 316-321,
            % 1999. https://doi.org/10.1115/1.2893982
            %

            if(~exist("gamma", "var"))
                gamma_m_R = 1.0;
                gamma_m_G = 1.0;
                gamma_J_R = 1.0;
                gamma_J_G = 1.0;
            else
                gamma_m_R = gamma(1);
                gamma_m_G = gamma(2);
                gamma_J_R = gamma(3);
                gamma_J_G = gamma(4);
            end
            
            m_R = obj.m_Rotor; % [kg], Rotor mass
            m_G = obj.m_Gen;   % [kg], Generator mass according to [?]
            
            J_R = obj.J_Rotor; % [kg-m^2], Rotor inertia according to [3]
            J_G = obj.J_Gen;   % [kg-m^2], Generator inertia according to [4]
            
            J_R = J_R*gamma_J_R;            J_G = J_G*gamma_J_G;
            m_R = m_R*gamma_m_R;            m_G = m_G*gamma_m_G;
            
            M       = zeros(42, 42);
            K_b     = zeros(42, 42);
            K_m     = zeros(42, 42);
            K_Omega = zeros(42, 42);
            G       = zeros(42, 42);
            
            M(1, 1) = m_R;              M(end - 2, end - 2) = m_G;
            M(2, 2) = m_R;              M(end - 1, end - 1) = m_G;
            M(3, 3) = J_R;              M(end    , end    ) = J_G;
            
            M_tmp = obj.input_shaft.inertia_matrix("full");
            K_tmp = obj.input_shaft.stiffness_matrix("full");
            
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
    end
    
end
