classdef Drivetrain
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
        stage         (1, :) Gear_Set;                                                                  % [-],      gearbox stages
        P_rated       (1, :)          {mustBeNumeric, mustBeFinite, mustBePositive} = 5.0e3;            % [kW],     Rated power
        n_rotor       (1, :)          {mustBeNumeric, mustBeFinite, mustBePositive} = 12.1;             % [1/min.], Rated rotor speed
        main_shaft    (1, :) Shaft                                                  = Shaft;            % [-],      Input Shaft
        m_Rotor       (1, :)          {mustBeNumeric, mustBeFinite, mustBePositive} = 110.0e3;          % [kg],     Rotor mass according to [3]
        J_Rotor       (1, :)          {mustBeNumeric, mustBeFinite, mustBePositive} = 57231535.0;       % [kg-m^2], Rotor mass moment of inertia according to [6]
        m_Gen         (1, :)          {mustBeNumeric, mustBeFinite, mustBePositive} = 1900.0;           % [kg],     Generator mass according to [4]
        J_Gen         (1, :)          {mustBeNumeric, mustBeFinite, mustBePositive} = 534.116;          % [kg-m^2], Generator mass moment of inertia [4]
        N_stg         (1, 1)          {mustBeNumeric, mustBeFinite, mustBePositive} = 3;                % [-],      Number of stages
        dynamic_model (1, :) string {mustBeMember(dynamic_model, ["Thomson_ToV", ...
                                                                  "Kahraman_1994", ...
                                                                  "Lin_Parker_1999"])} = "Thomson_ToV"; % which dynamic model should be used to perform modal analysis on the Drivetrain.
        gamma                scaling_factor;                                                            % Scaling factors
    end
    
    properties(Dependent)
        T_out;   % [N-m],    Output torque for each stage
        n_out;   % [1/min.], Output speed  for each stage
        u;       % [-],      Cumulative gear ratio
        S_H;     % [-],      Safety factor for surface durability (against pitting)
        S_F;     % [-],      Safety factor for bending strength
        S_shaft; % [-],      Safey factor for the shafts
    end
    
    properties(SetAccess = private)
        % to store the values of some dependent variables:
%         S_H_val     (1, :)          {mustBeNumeric, mustBeFinite, mustBePositive} = 1.25;   % [-],      Safety factor for surface durability (against pitting)
        S_H_val     (1, :)          {mustBeNumeric} = 1.25;   % [-],      Safety factor for surface durability (against pitting)
        S_F_val     (1, :)          {mustBeNumeric, mustBeFinite, mustBePositive} = 1.56;   % [-],      Safety factor for bending strength
        S_shaft_val (1, :)          {mustBeNumeric, mustBeFinite, mustBePositive} = 1.0;    % [-],      Safey factor for the shafts
    end
    
    methods
        function obj = Drivetrain(N_st, stage, P_r, n_r, inp_shaft, m_R, J_R, m_G, J_G)
            if(nargin == 0)
                N_st = 3;
                for idx = 1:N_st
                    stage(idx) = NREL_5MW.gear_stage(idx);
                end
                
                P_r = 5.0e3; % [kW], Rated power
                n_r = 12.1; % [1/min.], Input speed
                inp_shaft = Shaft;
                
                m_R = 110.0e3;      J_R = 57231535.0;
                m_G = 0.02*m_R;     J_G = 534.116;
                
            end
            
            obj.stage = stage;
            obj.P_rated = P_r;
            obj.n_rotor = n_r;
            obj.main_shaft = inp_shaft;
            
            obj.m_Rotor = m_R;          obj.J_Rotor = J_R;
            obj.m_Gen = m_G;            obj.J_Gen = J_G;
            
            obj.N_stg = N_st;
            
            try
                actxserver('KISSsoftCOM.KISSsoft');
                [obj.S_H_val, obj.S_F_val]     = obj.safety_factors_KS;
                [~          , obj.S_shaft_val] = obj.safety_factors;
            catch err
                warning(err.identifier, "%s", err.message)
                [obj.S_H_val, obj.S_shaft_val] = obj.safety_factors;
            end
            
            %      Gear stage               , Output shaft
            %      Normal module, Face width, Length, Diameter
            par = ["m_n1"       , "b_1"     , "d_1" , "L_1", ... % Stage 01
                   "m_n2"       , "b_2"     , "d_2" , "L_2", ... % Stage 02
                   "m_n3"       , "b_3"     , "d_3" , "L_3", ... % Stage 03
                   "J_R",                                    ... % M. M. Inertia (rotor)
                   "J_G",                                    ... % M. M. Inertia (generator)
                                                "d_s" , "L_s"]'; % Main shaft

            obj.gamma = scaling_factor(par, ones(size(par)));
        end
        
        function tab = disp(obj)
            %DISP display some properties of a Drivetrain object
            % description, symbol, unit, value
            if(isempty(obj))
                disp("\t0x0 empty Drivetrain object")
            else
                tab_str = ["Rated power",                                  "P",       "kW";
                    "Output Speed (Sun/Pinion)",                    "n_out",   "1/min.";
                    "Output Torque (Sun/Pinion)",                   "T_out",   "N-m";
                    "Minimum safety factor against pitting",        "S_Hmin",  "-";
                    "Safety factor against pitting (Sun/Pinion)",   "S_H1",    "-";
                    "Safety factor against pitting (Planet/Wheel)", "S_H2",    "-";
                    "Safety factor (Shaft)",                        "S",       "-";
                    "Type",                                         "-",       "-";
                    "Gear ratio",                                   "u",       "-";
                    "Number of planets",                            "p",       "-";
                    "Normal module",                                "m_n",     "mm";
                    "Normal pressure angle",                        "alpha_n", "deg.";
                    "Helix angle",                                  "beta",    "deg.";
                    "Face width",                                   "b",       "mm";
                    "Center distance",                              "a_w",     "mm";
                    "Number of teeth (Sun/Pinion)",                 "z_1",     "-";
                    "Number of teeth (Planet/Wheel)",               "z_2",     "-";
                    "Number of teeth (Ring)",                       "z_3",     "-";
                    "Profile shift coefficient (Sun/Pinion)",       "x_1",     "-";
                    "Profile shift coefficient (Planet/Wheel)",     "x_2",     "-";
                    "Profile shift coefficient (Ring)",             "x_3",     "-";
                    "Reference diameter (Sun/Pinion)",              "d_1",     "mm";
                    "Reference diameter (Planet/Wheel)",            "d_2",     "mm";
                    "Reference diameter (Ring)",                    "d_3",     "mm";
                    "Mass (Sun/Pinion)",                            "m_1",     "kg";
                    "Mass (Planet/Wheel)",                          "m_2",     "kg";
                    "Mass (Ring)",                                  "m_3",     "kg";
                    "Mass mom. inertia (Sun/Pinion)",               "J_xx1",   "kg-m^2";
                    "Mass mom. inertia (Planet/Wheel)",             "J_xx2",   "kg-m^2";
                    "Mass mom. inertia (Ring)",                     "J_xx3",   "kg-m^2";
                    "Mass mom. inertia (Sun/Pinion)",               "J_yy1",   "kg-m^2";
                    "Mass mom. inertia (Planet/Wheel)",             "J_yy2",   "kg-m^2";
                    "Mass mom. inertia (Ring)",                     "J_yy3",   "kg-m^2";
                    "Mass mom. inertia (Sun/Pinion)",               "J_zz1",   "kg-m^2";
                    "Mass mom. inertia (Planet/Wheel)",             "J_zz2",   "kg-m^2";
                    "Mass mom. inertia (Ring)",                     "J_zz3",   "kg-m^2";
                    ];
                
                Parameter = tab_str(:, 1);
                Symbol    = tab_str(:, 2);
                Unit      = tab_str(:, 3);
                
                tab_val = cell(obj.N_stg + 2, 1);
                
                tab_val{1} = table(Parameter, Symbol);
                
                for idx = 1:obj.N_stg
                    val_stg = {obj.P_rated;
                        obj.n_out(idx);
                        obj.T_out(idx);
                        1.25;
                        obj.S_H(2*idx - 1);
                        obj.S_H(2*idx);
                        obj.S_shaft(idx + 1);
                        obj.stage(idx).configuration;
                        obj.stage(idx).u;
                        obj.stage(idx).N_p;
                        obj.stage(idx).m_n;
                        obj.stage(idx).alpha_n;
                        obj.stage(idx).beta;
                        obj.stage(idx).b;
                        obj.stage(idx).a_w;
                        obj.stage(idx).z(1);
                        obj.stage(idx).z(2);
                        "-*-";
                        obj.stage(idx).x(1);
                        obj.stage(idx).x(2);
                        "-*-";
                        obj.stage(idx).d(1);
                        obj.stage(idx).d(2);
                        "-*-";
                        obj.stage(idx).mass(1);
                        obj.stage(idx).mass(2);
                        "-*-";
                        obj.stage(idx).J_x(1);
                        obj.stage(idx).J_x(2);
                        "-*-";
                        obj.stage(idx).J_y(1);
                        obj.stage(idx).J_y(2);
                        "-*-";
                        obj.stage(idx).J_z(1);
                        obj.stage(idx).J_z(2);
                        "-*-";
                        };
                    
                    if(strcmp(obj.stage(idx).configuration, "planetary"))
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
                    fprintf("Gear stages:\n");
                    disp(tab);
                    fprintf("Main shaft:\n");
                    disp(obj.main_shaft);
                    clear tab;
                end
            end            
        end
        
        function tab = comparison(ref, sca)
            tab_stg = cell(ref.N_stg, 1);
            
            for idx = 1:(ref.N_stg)
                [~, tab_str] = stage_comparison(ref, sca, idx);
                tab_tmp = table(tab_str(:, 4), tab_str(:, 5), tab_str(:, 6), 'variableNames', ["Reference", "Scale", "Ratio"]);
                tab_stg{idx} = table(tab_tmp, 'variableNames', sprintf("Stage_%d", idx));
            end
            
            tab_left =  table(tab_str(:, 1), tab_str(:, 2), 'variableNames', ["Parameter", "Symbol"]);
            tab_right = table(tab_str(:, 3), 'variableNames', "Unit");
            
            tab = [tab_left, tab_stg{:}, tab_right];
            
            if(nargout == 0)
                disp(tab);
                clear tab;
            end
            
        end
        
        function [tab, tab_str] = stage_comparison(ref, sca, idx)
            if(ref.N_stg ~= sca.N_stg)
                error("different number of gear stages.");
            elseif(idx > ref.N_stg)
                error("idx is bigger than the number of stages.");
            elseif(idx < 1)
                error("idx is smaller than 1.");
            end
            
            tab_str = {"Rated power",                                  "P",       "kW",     ref.P_rated,                sca.P_rated,                sca.P_rated               /ref.P_rated;                % 1
                       "Output Speed (Sun/Pinion)",                    "n_out",   "1/min.", ref.n_out(idx),             sca.n_out(idx),             sca.n_out(idx)            /ref.n_out(idx);             % 2
                       "Output Torque (Sun/Pinion)",                   "T_out",   "N-m",    ref.T_out(idx),             sca.T_out(idx),             sca.T_out(idx)            /ref.T_out(idx);             % 3
                       "Safety factor against pitting (Sun/Pinion)",   "S_H1",    "-",      ref.S_H(2*idx - 1),         sca.S_H(2*idx - 1),         sca.S_H(2*idx - 1)        /ref.S_H(2*idx - 1);         % 4
                       "Safety factor against pitting (Planet/Wheel)", "S_H2",    "-",      ref.S_H(2*idx),             sca.S_H(2*idx),             sca.S_H(2*idx)            /ref.S_H(2*idx);             % 5
                       "Safety factor (Shaft)",                        "S",       "-",      ref.S_shaft(idx + 1),       sca.S_shaft(idx + 1),       sca.S_shaft(idx + 1)      /ref.S_shaft(idx + 1);
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
            for idx = 1:obj.N_stg
                subplot(1, obj.N_stg, idx)
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
            color = linspecer(6, "qualitative");
            
            hold on;
            C_s = [obj.main_shaft.L/2.0 0.0]' + C_0;
            rectangle(obj.main_shaft, C_s, color(6, :), "edgeColor", "k", "lineStyle", "-" , "faceColor", color(6, :));
            
            C_x = [obj.main_shaft.L + C_0(1);
                   obj.stage(1).carrier.b + obj.stage(1).out_shaft.L;
                   obj.stage(2).carrier.b + obj.stage(2).out_shaft.L];

            for idx = 1:obj.N_stg
                rectangle(obj.stage(idx), [sum(C_x(1:idx)) C_0(2)]');
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
            
            flag_im = any(imag(f_n) ~= 0.0);
            if(flag_im)
                idx = (imag(f_n) ~= 0.0);
                f_n(idx) = 0.0;
                
                f_n = [f_n(~idx);
                       f_n(idx)];
                   
                mode_shape = [mode_shape(:, ~idx), mode_shape(:, idx)];
            end
            
            flag_RB = any(abs(f_n) < 1.0e-2);
            if(flag_RB)
                idx = (abs(f_n) < 1.0e-2);
                f_n(idx) = 0.0;
                
                f_n = [f_n(~idx);
                       f_n(idx)];
                   
                mode_shape = [mode_shape(:, ~idx), mode_shape(:, idx)];
            end
            
            % Normalizing the mode shapes so that the maximum is always +1:
            for idx = 1:length(f_n)
                [ms_max, n] = max(abs(mode_shape(:, idx)));
                mode_shape(:, idx) = mode_shape(:, idx)*sign(mode_shape(n, idx))./ms_max;
            end
            
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
            
            N_DOF = ones(obj.N_stg + 1, 1)*2.0;
            
            for idx = 1:obj.N_stg
                if(strcmp(obj.stage(idx).configuration, "parallel"))
                    tmp = obj.stage(idx).N_p + 1;
                elseif(strcmp(obj.stage(idx).configuration, "planetary"))
                    tmp = obj.stage(idx).N_p + 2;
                end
                
                N_DOF(idx + 1) = N_DOF(idx) + tmp;
            end
            
            n = N_DOF(1:end - 1);
            N_DOF = N_DOF(end);
            
            M = zeros(N_DOF, N_DOF);
            K = zeros(N_DOF, N_DOF);
            
            M(1, 1) = J_R;              M(end, end) = J_G;
            
            M_tmp = obj.main_shaft.inertia_matrix("torsional");
            K_tmp = obj.main_shaft.stiffness_matrix("torsional");
            
            M(1:2, 1:2) = M(1:2, 1:2) + M_tmp;        K(1:2, 1:2) = K(1:2, 1:2) + K_tmp;
            
            for idx = 1:obj.N_stg
                [M_tmp, K_tmp] = obj.stage(idx).Kahraman_1994;
                
                jdx = n(idx);
                kdx = jdx:(jdx + length(M_tmp) - 1);
                
                M(kdx, kdx) = M(kdx, kdx) + M_tmp;
                K(kdx, kdx) = K(kdx, kdx) + K_tmp;
            end
            
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
            
            N_DOF = 6;
            
            for idx = 1:obj.N_stg
                if(strcmp(obj.stage(idx).configuration, "parallel"))
                    tmp = obj.stage(idx).N_p + 1;
                elseif(strcmp(obj.stage(idx).configuration, "planetary"))
                    tmp = obj.stage(idx).N_p + 2;
                end
                
                N_DOF = N_DOF + tmp*3;
            end
            
            M       = zeros(N_DOF, N_DOF);
            K_b     = zeros(N_DOF, N_DOF);
            K_m     = zeros(N_DOF, N_DOF);
            K_Omega = zeros(N_DOF, N_DOF);
            G       = zeros(N_DOF, N_DOF);
            
            M(1, 1) = m_R;              M(end - 2, end - 2) = m_G;
            M(2, 2) = m_R;              M(end - 1, end - 1) = m_G;
            M(3, 3) = J_R;              M(end    , end    ) = J_G;
            
            M_tmp = obj.main_shaft.inertia_matrix("full");
            K_tmp = obj.main_shaft.stiffness_matrix("full");
            
            idx = [1 5 6 7 11 12];
            M_tmp(idx, :) = [];                 M_tmp(:, idx) = [];
            K_tmp(idx, :) = [];                 K_tmp(:, idx) = [];
            
            M(  1:6, 1:6) = M(  1:6, 1:6) + M_tmp;
            K_b(1:6, 1:6) = K_b(1:6, 1:6) + K_tmp;
            
            for idx = 1:obj.N_stg
                [M_tmp, ~, K_b_tmp, K_m_tmp, K_Omega_tmp, G_tmp] = obj.stage(idx).Lin_Parker_1999;
                
                jdx = 15*idx - 11;
                if(strcmp(obj.stage(idx).configuration, "parallel"))
                    kdx = jdx:(jdx + 8);
                elseif(strcmp(obj.stage(idx).configuration, "planetary"))
                    kdx = jdx:(jdx + 17);
                end
                
                M(kdx, kdx) = M(kdx, kdx) + M_tmp;
                
                K_b(    kdx, kdx) = K_b(    kdx, kdx) + K_b_tmp;
                K_m(    kdx, kdx) = K_m(    kdx, kdx) + K_m_tmp;
                K_Omega(kdx, kdx) = K_Omega(kdx, kdx) + K_Omega_tmp;
                
                G(kdx, kdx) = G(kdx, kdx) + G_tmp;
                
            end
            
            K = @(Om)(K_b + K_m - K_Omega*Om^2);
        end
        
        function SIMPACK_version(obj, varargin)
            if(nargin == 1)
                version = "2018";
                task = "lsa";
            elseif(nargin == 2)
                version = "2018";
                task = varargin{1};
            elseif(nargin == 3)
                version = varargin{1};
                task = varargin{2};
            else
                error("too many arguments.");
            end
            
            % create SIMPACK COM solver interface:
            solver = actxserver(sprintf('Simpack.Slv.%s', version));
            
            % Open model:
            file = dir(sprintf("@%s\\*.spck", class(obj)));
            file_name = sprintf("%s\\%s", file.folder, file.name);
            model = solver.Spck.openModel(file_name);
            
            switch(task)
                case "lsa"
                    % Linear System Analysis:
                    solver.Spck.Slv.lsa(model);
                    
                    % the following part seems useless:
                    % create SIMPACK COM result interface:
                    post = actxserver(sprintf("Simpack.Post.%s", version));

                    % add project:
                    project = post.Spck.addProject();
                    file = dir(sprintf("@%s\\*\\*.lsa.sbr", class(obj)));
                    file_name = sprintf("%s\\%s", file.folder, file.name);

                    % add Linear System Response result file to project:
                    result = project.addResultFile(file_name);

                case "eigen"
                    % Modal Analysis:
                    result = solver.Spck.Slv.eigen(model, false);
                    nx = result.numEigenvalues; % number of eigenvalues
                    freq = zeros(nx, 1);

                    for idx = 1:nx
                        freq(idx) = result.freq(idx - 1);
                    end

                    idx = find(freq > 1.0e-5, 1, 'first');
                    freq = freq(idx:end);
                    
                case "ssm"
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
        
        function [SH_vec, SF_vec] = safety_factors_KS(obj, varargin)
            
            if(nargin == 1)
                save_report = false;
                show_report = false;
            else
                save_report = varargin{1}; % true
                show_report = varargin{2}; % true
            end            
            
            S_Hmin = 1.25;      % [-],  Minimum required safety factor for surface durability according to IEC 61400-4.
            S_Fmin = 1.56;      % [-],  Minimum required safety factor for surface durability according to IEC 61400-4.
            L_h    = 20*365*24; % [h],  Required life
            K_A    = 1.25;      % [-],  Application factor
            
            max_Np = max([obj.stage.N_p]) + 1;
            
            SH_vec = zeros(obj.N_stg, max_Np);
            SF_vec = zeros(obj.N_stg, max_Np);
            
            speed = [obj.n_rotor; obj.n_out];
            
            for idx = 1:obj.N_stg
                ks = obj.stage(idx).toKISSsoft();
                
                ks.SetVar('ZS.P'       , num2str(obj.P_rated, '%.6f'));    % rated power
                ks.SetVar('ZS.Pnominal', num2str(obj.P_rated, '%.6f'));    % rated power
                
                flag_impul = true;
                if(strcmp(obj.stage(idx).configuration, "parallel"))
                    flag_impul = true;
                    template_code = 12;
                elseif(strcmp(obj.stage(idx).configuration, 'planetary'))
                    n_c = speed(idx);
                    T_c =obj.P_rated/(n_c*pi/30.0);
                    n_planet = n_c*(1.0 - obj.stage(idx).z(3)./obj.stage(idx).z(2));
                    ks.SetVar('ZR[1].n'               , num2str(abs(n_planet), '%.6f')); % speed at gear 2: planet
                    ks.SetVar('ZR[2].n'               , num2str(0.0          , '%.6f')); % speed at gear 3: ring
                    ks.SetVar('ZS.Planet.nStegnominal', num2str(n_c          , '%.6f')); % speed at planet carrier
                    ks.SetVar('ZS.Planet.nSteg'       , num2str(n_c          , '%.6f')); % speed at planet carrier
                    
                    ks.SetVar('ZR[1].T'               , num2str(abs(n_planet), '%.6f')); % speed at gear 2: planet
                    ks.SetVar('ZS.Planet.TStegnominal', num2str(T_c          , '%.6f')); % torque at planet carrier
                    ks.SetVar('ZS.Planet.nSteg'       , num2str(T_c          , '%.6f')); % torque at planet carrier
                    
                    flag_impul = false;
                    template_code = 14;
                end
                
                n_1 = obj.n_out(idx);
                T_1 = obj.P_rated/(n_1*pi/30.0);
                ks.SetVar('ZR[0].Tnominal', num2str(T_1, '%.6f')); % torque at gear 1
                ks.SetVar('ZR[0].n'       , num2str(n_1, '%.6f')); % speed at gear 1
                ks.SetVar('ZP[0].Impuls'  , num2str(flag_impul));  % sense of rotation, gear 1
                
                ks.SetVar('ZS.SSi.Flanke'      , num2str(S_Hmin, '%.6f'));
                ks.SetVar('ZS.SSi_0to05.Flanke', num2str(S_Hmin, '%.6f'));
                ks.SetVar('ZS.SSi_at1.Flanke'  , num2str(S_Hmin, '%.6f'));
                ks.SetVar('ZS.SSi_2to20.Flanke', num2str(S_Hmin, '%.6f'));
                ks.SetVar('ZS.SSi_fix.Flanke'  , num2str(S_Hmin, '%.6f'));
                
                ks.SetVar('ZS.SSi.Fuss'      , num2str(S_Fmin, '%.6f'));
                ks.SetVar('ZS.SSi_0to05.Fuss', num2str(S_Fmin, '%.6f'));
                ks.SetVar('ZS.SSi_at1.Fuss'  , num2str(S_Fmin, '%.6f'));
                ks.SetVar('ZS.SSi_2to20.Fuss', num2str(S_Fmin, '%.6f'));
                ks.SetVar('ZS.SSi_fix.Fuss'  , num2str(S_Fmin, '%.6f'));
                
                ks.SetVar('ZS.H'         , num2str(L_h, '%.6f'));
                ks.SetVar('ZS.KA'        , num2str(K_A, '%.6f'));
                
                n_gear = numel(obj.stage(idx).z);
                
                % [5], Table 3, p. 35:
                switch(obj.stage(idx).N_p)
                    case 3
                        K_gam = 1.1;
                    case 4
                        K_gam = 1.25;
                    case 5
                        K_gam = 1.35;
                    case 6
                        K_gam = 1.44;
                    case 7
                        K_gam = 1.47;
                    otherwise
                        K_gam = 1.0;
                end
                
                ks.SetVar('ZS.Kgam', num2str(K_gam, '%.6f'));
                
                if(~ks.CalculateRetVal())
                    ks.ReleaseModule();
                    error("Error in KISSsoft calculation.");
                end
                
                kv = str2double(ks.GetVar('ZS.KVcalc'));
                flag_kv = false;
                if(kv < 1.05)
                    flag_kv = true;
                    for jdx = 1:n_gear
                        ks.SetVar(sprintf('ZP[%d].KV.KV', jdx), num2str(1.05, '%.6f'));
                    end
                    
                    ks.SetVar('Zst.KVFlag', num2str(flag_kv));
                end
                
                flag_khb = false;
                for jdx = 1:2
                    khb = str2double(ks.GetVar(sprintf('ZP[%d].KHb', jdx - 1)));
                    
                    if(khb < 1.15)
                        flag_khb = true;
                        ks.SetVar(sprintf('ZP[%d].KHb', jdx - 1), num2str(1.15, '%.6f'));
                    end
                end
                
                if(flag_khb)
                    ks.SetVar('Zst.KHbVariant', num2str(true));
                end
                
                if(~ks.CalculateRetVal())
                    ks.ReleaseModule();
                    error("Error in KISSsoft calculation.");
                end
                
                if(flag_kv || flag_khb)
                    ks.Calculate();
                end

                for jdx = 1:n_gear
                    SH_vec(jdx, idx) = str2double(ks.GetVar(sprintf('ZPP[%d].Flanke.SH', jdx - 1)));
                    SF_vec(jdx, idx) = str2double(ks.GetVar(sprintf('ZPP[%d].Fuss.SF'  , jdx - 1)));
                end
                
                if(save_report)
                    if(strcmp(class(obj), "Drivetrain"))
                        report_name = sprintf("%s\\stage_%02d.rtf" , pwd, idx);
                        file_name   = sprintf("%s\\stage_%02d.Z0%d", pwd, idx, template_code);
                    else
                        report_name = sprintf("%s\\@%s\\stage_%02d.rtf" , pwd, class(obj), idx);
                        file_name   = sprintf("%s\\@%s\\stage_%02d.Z0%d", pwd, class(obj), idx, template_code);
                    end
                    
                    ks.SaveFile(file_name);

                    ks.ReportWithParameters(sprintf("C:\\Program Files (x86)\\KISSsoft 03-2017\\rpt\\Z0%dLe0.RPT", template_code), ... % template file
                                            report_name, show_report, 0); % output format
                else
                    ks.Report(show_report);
                end
                                    
                ks.ReleaseModule();
                clear("ks");
            end
            
            SH_vec(SH_vec == 0) = [];
            SF_vec(SF_vec == 0) = [];
            
            if(isrow(SH_vec))
                SH_vec = SH_vec';
            end
            
            if(isrow(SF_vec))
                SF_vec = SF_vec';
            end
        end
        
        %% Scaling:
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
            
            if((length(gamma) ~= length(obj_ref.gamma)) || (~isa(gamma, "scaling_factor")))
                error("gamma must be a scaling_factor object with %d elements.", length(obj_ref.gamma));
            end
            
            key_set = obj_ref.gamma.name;
            
            gamma_mR  = 1.0; %gamma("m_R");
            gamma_JR  = gamma("J_R");
            gamma_mG  = 1.0; %gamma("m_G");
            gamma_JG  = gamma("J_G");
%             gamma_Ls  = gamma("L_s");
%             gamma_ds  = gamma("d_s");
            
            gamma_Torque = gamma_P/gamma_n; % Applied torque
            
            % Main shaft:
%             gamma_Ls = gamma("L_3"); % length
            gamma_ds = nthroot(gamma_Torque, 3.0);
            
            % Shaft diameters are scaled in the same way
            idx_dd = contains(key_set, "d");
            idx_LL = contains(key_set, "L");
            
            gamma(key_set(idx_dd)) = gamma_ds*ones(sum(idx_dd), 1);
            
            gamma_dd = gamma(key_set(idx_dd));
            gamma_LL = gamma(key_set(idx_LL));
            
            stage_sca = [Gear_Set, Gear_Set, Gear_Set];
            
            for idx = 1:obj_ref.N_stg
                jdx = 4*idx + (-3:0);
                gamma_stage = gamma(key_set(jdx));
                
                stage_sca(idx) = obj_ref.stage(idx).scale_aspect(gamma_stage, "Gear_Set");
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
        
        function obj_sca = scale_aspect(obj_ref, gamma_P, gamma_n, gm_val, aspect)
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
                error("Aspect contains NO parameter related to %s.", class(obj_ref));
            end
            
            lin = size(aspect, 1);
            
            if(numel(gm_val) ~= lin)
                error("gm_val should have the same number of lines of aspect.")
            end
            
            for idx = 1:lin
                asp_tmp = aspect(idx, :);
                asp_tmp = asp_tmp(asp_tmp ~= "*");
                col = numel(asp_tmp);
                
                for jdx = 1:col
                    gamma_full(aspect(idx, jdx)) = gm_val(idx, :);
                end
            end
            
            % Scaling the shaft's diameter:
            gamma_T = gamma_P/gamma_n; % Applied torque scaling factor
            gamma_ds = nthroot(gamma_T, 3.0);

            idx_D = contains(key_set, "d");
            gamma_full(key_set(idx_D)) = gamma_ds*ones(sum(idx_D), 1);
            
            if(strcmp(class(obj_ref), "Drivetrain"))
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
            % ["L_1", "L_2", "L_3", "L_s"; ...
            %  "J_R", "J_G", "*",   "*" ]; % '*' are used to complete the
            %  string array.
            %
            
            % [TODO]:
            % "K_MMI_stage": scales parameters related to the stiffness and
            % mass moment of inertia, including the gear stages:
            % min |1 - f(sca)/f(ref)|^2 = f(x)
            % subjected to: S(min) - S(sca) <= 0 = g(x)
            %               S(ref) - S(sca)  = 0 = h(x)
            % ["m_n1", "b_1"; ...
            %  "m_n2", "b_2"; ...
            %  "m_n3", "b_3"; ...
            %  "L_1" , "L_2", "L_3"; ...
            %  "J_R", "J_G"; ...
            %  "L_s"];
            % "K_MMI_gear"
            % "M_MMI_stage_det"
            % "K_MMI_gear_det"
            % "K_shaft" : scales only parameters related to the shaft's 
            % "K": scales only stiffness-related parameters. The mesh
            % stiffnesss is assumed to be proportinal to the gear's face
            % width.

        end
        
        function [obj_sca, gamma_val, res, gamma_sep] = scaled_version(obj_ref, P_scale, n_R_scale, normalize_freq, N_freq, aspect_set, varargin)
            %SCALED_VERSION returns an scaled object together with the 
            % scaling factor for its parameters and residual error from the
            % scaling optimization process. Different scaling can be
            % obtained by using:
            % - aspect_set: Different model approaches;
            % - normalize_freq: Normalizing or not the target resonances;
            % - N_freq: choosing how many resonances should be considered
            % during the scaling optimization process;
            % Additionaly, one can define an initial value for gamma to be
            % used on the optimization process
            
            if(~isempty(varargin))
                gamma_prev_val = varargin{1};
            else
                gamma_prev_val = ones(size(obj_ref.gamma.name))*0.5;
            end
            
            gamma_prev = scaling_factor(obj_ref.gamma.name, gamma_prev_val);
            
            % Input scaling factors:
            gm_P = P_scale/obj_ref.P_rated;
            gm_n = n_R_scale/obj_ref.n_rotor;
            gm_T = gm_P/gm_n;
            
            % Scaling factors from dimensional analysis:
            gamma_length = nthroot(gm_T,     3.0);
            gamma_MMI    = power(  gm_P, 5.0/3.0);
            
            S_H_ref = obj_ref.S_H;
            f_n_ref = obj_ref.resonances(N_freq, normalize_freq);
            
            if(strcmp(class(obj_ref), "Drivetrain"))
                fprintf("Scaling drivetrain to rated power %.1f kW, which is %.2f %% of its reference.\n", gm_P*[obj_ref.P_rated 100.0]);
            else
                fprintf("Scaling %s drivetrain to rated power %.1f kW, which is %.2f %% of its reference.\n", class(obj_ref), gm_P*[obj_ref.P_rated 100.0]);
            end
            
            id_1 = "prog:input";
            id_2 = "MATLAB:nearlySingularMatrix";
            id_3 = "MATLAB:COM:InvalidProgid";
            warning("off", id_1);
            warning("off", id_2);
            warning("off", id_3);
            
            %% 1. Stage scaling:
            aspect_1 = "stage";
            
            if(any(fields(aspect_set) == aspect_1))
                fprintf("Optimizing w.r.t. [%s]...\n", upper(aspect_1));
                
                gm_01 = mean(gamma_prev(["m_n1" "b_1" "m_n2" "b_2" "m_n3" "b_3"]));
                gm_val_1  = zeros(obj_ref.N_stg, 1);
                res_1 = zeros(obj_ref.N_stg, 1);

                for idx = 1:obj_ref.N_stg
                    jdx = 2*idx + (-1:0);
                    n_stage = obj_ref.n_out(idx)*gm_n;

                    [~, gm_val_1(idx, :), res_1(idx, :)] = obj_ref.stage(idx).scaled_version(P_scale, n_stage, S_H_ref(jdx), aspect_1, gm_01);
                    gm_01 = gm_val_1(idx, :);
                end
            elseif(any(fields(aspect_set) == "DA"))
                m_n_tmp =[obj_ref.stage.m_n]';
                gm_val_1 = Rack.module(m_n_tmp*gamma_length, "calc", "nearest")./m_n_tmp;
                
                res_1 = Inf(3, 1);
            else
                gm_val_1 = [mean(gamma_prev(["m_n1" "b_1"]));
                            mean(gamma_prev(["m_n2" "b_2"]));
                            mean(gamma_prev(["m_n3" "b_3"]))];
                       
                res_1 = Inf(3, 1);
            end
            
            gm_val_1  =  gm_val_1([1 1 2 2 3 3]);
            res_1     = res_1([1 1 2 2 3 3]);
            
            %% 2. Gears:
            aspect_2 = "gear";
            
            if(any(fields(aspect_set) == aspect_2))
                fprintf("Optimizing w.r.t. [%s]...\n", upper(aspect_2));
                
                gm_02     = gm_val_1([1 1]);
                gm_val_2  = zeros(2*obj_ref.N_stg, 1);
                res_2     = zeros(2*obj_ref.N_stg, 1);

                for idx = 1:obj_ref.N_stg
                    jdx = 2*idx + (-1:0);
                    n_stage = obj_ref.n_out(idx)*gm_n;

                    [~, gm_val_2(jdx, :), res_2(jdx, :)] = obj_ref.stage(idx).scaled_version(P_scale, n_stage, S_H_ref(jdx), aspect_2, gm_02);
                    gm_02 = gm_val_2(jdx, :);
                end
            elseif(any(fields(aspect_set) == "DA"))
                m_n_tmp = [obj_ref.stage.m_n]';
                gamma_mn = Rack.module(m_n_tmp*gamma_length, "calc", "nearest")./m_n_tmp;
                
                gm_val_2          = ones(6, 1)*gamma_length;
                gm_val_2(1:2:end) = gamma_mn;
                
                res_2 = Inf(6, 1);
            else
                gm_val_2 = gamma_prev(["m_n1" "b_1" "m_n2" "b_2" "m_n3" "b_3"]);
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
            
            gamma_sca = obj_ref.gamma;
            
            gamma_sca(["m_n1" "b_1" ...
                       "m_n2" "b_2" ...
                       "m_n3" "b_3"]) = gm_val_12;
                 
            if(strcmp(class(obj_ref), "Drivetrain"))
                obj_12 = scale_all(gm_P, gm_n, gamma_sca);
            else
                scale_func = str2func(class(obj_ref));
                obj_12 = scale_func(gm_P, gm_n, gamma_sca);
            end
            
            %% 3. Shaft stiffness and Mass moment of inertia of rotor and generator:
            aspect_3 = "K_MMI";
            
            opt_solver = optimoptions("fmincon", "display", "notify");

            if(any(fields(aspect_set) == aspect_3))
                fprintf("Optimizing w.r.t. [%s]...\n", upper(aspect_3));
                
                asp_set = aspect_set.(aspect_3);
                n = size(asp_set, 1);
                
                fun_asp  = @(x)(1.0 - obj_12.scale_resonances(N_freq, normalize_freq, gm_P, gm_n, x, asp_set)./f_n_ref);
                fun_min = @(x)(norm(fun_asp(x))^2);

                gamma_min = ones(n, 1)*1.0e-6;
                gamma_Max = ones(n, 1);
                
                idx_L = contains(obj_ref.gamma.name, "L");
                idx_J = contains(obj_ref.gamma.name, "J");

                gm_03 = [mean(gamma_prev(idx_L)); ... % length
                         mean(gamma_prev(idx_J))]; % mass mom. inertia

                constraint_fun = @(x)deal([], fun_asp(x)); % inequalities, equalities

                [gm_val_3, res_3, ~] = fmincon(fun_min, gm_03, [], [], [], [], gamma_min, gamma_Max, constraint_fun, opt_solver);
            elseif(any(aspect_set == "DA"))
                gm_val_3 = [gamma_length;
                            gamma_MMI];
                res_3 = Inf;
            else
                gm_val_3 = [mean(gamma_prev(["L_1" "L_2" "L_3" "L_s"]));
                            mean(gamma_prev(["J_R" "J_G"]))];
                       
                res_3 = Inf;
            end
            
            gamma_sca(["L_1" "L_2" "L_3" "L_s"]) = gm_val_3(1)*ones(4, 1);
            gamma_sca(["J_R" "J_G"])             = gm_val_3(2)*ones(2, 1);

            %% 4. Detailed version of 3:
            aspect_4 = "K_MMI_det";
            
            if(any(fields(aspect_set) == aspect_4))
                fprintf("Optimizing w.r.t. [%s]...\n", upper(aspect_4));
                
                asp_set = aspect_set.(aspect_4);
                n = size(asp_set, 1);
                
                % function for aspect_4: gamma_P and gamma_n are set to 1.0
                % because the drivetrain is already scaled for these
                % parameters.
                fun_asp  = @(x)(1.0 - obj_12.scale_resonances(N_freq, normalize_freq, gm_P, gm_n, x, asp_set)./f_n_ref);
                fun_min = @(x)(norm(fun_asp(x))^2);

                gamma_min = ones(n, 1)*1.0e-6;
                gamma_Max = ones(n, 1);

                gm_04 = ones(n, 1);
                
                for idx = 1:n
                    gm_04(idx) = mean(gamma_sca(asp_set(idx, :)));
                end

                constraint_fun = @(x)deal([], fun_asp(x)); % inequalities, equalities

                [gm_val_4, res_4, ~] = fmincon(fun_min, gm_04, [], [], [], [], gamma_min, gamma_Max, constraint_fun, opt_solver);
            elseif(any(fields(aspect_set) == "DA"))
                gm_val_4      = gamma_length*ones(6, 1);
                gm_val_4(4:5) = gamma_MMI; 
                
                res_4 = Inf;
            else
                gm_val_4 = gamma_prev(["L_1" "L_2" "L_3" "L_s" "J_R" "J_G"]);
                res_4    = Inf;
            end
            
            if(all(isinf([res_3 res_4])) || (res_3 <= res_4))
                gm_val_34 = gm_val_3;
                res_34 = res_3;
                
                gamma_sca(["L_1" "L_2" "L_3" "L_s"]) = gm_val_34(1)*ones(4, 1);
                gamma_sca(["J_R" "J_G"            ]) = gm_val_34(2)*ones(2, 1);
            else
                gm_val_34 = gm_val_4;
                res_34 = res_4;
                
                for idx = 1:n
                    asp_tmp = asp_set(idx, :);
                    asp_tmp(asp_tmp == "*") = [];
                    
                    m = size(asp_tmp);
                    gamma_sca(asp_set(idx, :)) = gm_val_34(idx)*ones(m);
                end
            end
            
            %% 5. Drivetrain scaling:
            aspect_5 = "Drivetrain";
            
            if(any(fields(aspect_set) == aspect_5))
                fprintf("Optimizing w.r.t. [%s]...\n", upper(aspect_5));
                
                fprintf("\n\tto be done later...\n\n");
            else
                
            end
            
            %% Post processing:
            % Analysis of the residuals at each step:
            res_pp = [mean(res_1), mean(res_2),    res_3,    res_4];
            idx = ~isinf(res_pp);
            res_pp = res_pp(idx);
            
            asp_name = fields(aspect_set);
            
            if(numel(asp_name) ~= sum(idx))
                asp_name = asp_name(idx);
            end
            
            fprintf("Scale: %.1f kW = %.2f %% of Ref.\n", gm_P*[obj_ref.P_rated, ...
                                                                100.0]);

            if(isempty(res_pp))
                fprintf("\tScaled version obtained by dimensional analysis scaling rules.\n");
            else
                [~, sorted_idx] = sort(res_pp);
                
                fprintf("\tOptimization residua:\n");

                for idx = sorted_idx
                    fprintf("\t%d. %s:\t%.5e\n", idx, asp_name{idx}, res_pp(idx));
                end
            end
            
            gm_d = nthroot(gm_P/gm_n, 3.0);
            
            gamma_sca(["d_1" "d_2" "d_3" "d_s"]) = gm_d*ones(4, 1);
            gamma_val = gamma_sca.value;
            
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
            
            if(strcmp(class(obj_ref), "Drivetrain"))
                obj_sca = scale_all(gm_P, gm_n, gamma_sca);
            else
                scale_func = str2func(class(obj_ref));
                obj_sca = scale_func(gm_P, gm_n, gamma_sca);
            end
            
            warning("on", id_1);
            warning("on", id_2);
            warning("on", id_3);
        end
        
        function [gamma, res, SH, f_n, mode_shape, k_mesh, gamma_asp] = scaled_sweep(obj_ref, P_scale, n_R_scale, normalize_freq, N_freq, aspect_set)
            %SCALED_SWEEP performs a sweep on the rated power parameter of
            % the object. Returns the scaling factors gamma
            %
            % see also SCALED_VERSION
            %
            
            if(strcmp(class(obj_ref), "Drivetrain"))
                fprintf("Reference Drivetrain with rated power %.1f kW.\n", obj_ref.P_rated);
            else
                fprintf("Reference %s drivetrain with rated power %.1f kW.\n", class(obj_ref), obj_ref.P_rated);
            end
            
            n_P = numel(P_scale);
            n_fn = numel(obj_ref.modal_analysis);
            
            gamma = zeros(length(obj_ref.gamma), n_P);
            res = struct;
            gamma_asp = struct;
            SH = zeros(numel(obj_ref.S_H), n_P);
            f_n = zeros(n_fn, n_P);
            mode_shape = zeros(n_fn, n_fn, n_P);
            k_mesh = zeros(numel(obj_ref.stage), n_P);
            
            gm_P = P_scale./obj_ref.P_rated;
            
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
            
            tick_x = ["R"  , ...
                      "c_1", "p_{11}", "p_{12}", "p_{13}", "s_1", ...
                      "c_2", "p_{21}", "p_{22}", "p_{23}", "s_2", ...
                      "W_3", "P_3", ...
                      "G"];
            idx_tick = 2:2:14;
            
            font_setting  = {"fontName", "Times", "fontSize", 12.0};
            latex_setting = {"fontName", "Times", "fontSize", 12.0, "interpreter", "LaTeX"};
            
            fig_axes = findobj(gcf, "Type", "Axes");
            
            set(get(fig_axes(1), "ylabel"), "string", "$S_H$"             , latex_setting{:});
            set(get(fig_axes(5), "ylabel"), "string", "$\theta_i$, [deg.]", latex_setting{:});
            set(get(fig_axes(6), "ylabel"), "string", "$\theta_i$, [deg.]", latex_setting{:});

            set(fig_axes(2:6), 'xlim'      , [1 14]);
            set(fig_axes(2:6), 'ylim'      , [-1 1]*100);
            set(fig_axes(2:6), 'xtick'     ,        idx_tick);
            set(fig_axes(2:6), 'xticklabel', tick_x(idx_tick));
            set(fig_axes     , font_setting{:});
            
            label_figure(fig_axes, font_setting);
            
            fig_name = @(i)(sprintf("plots\\sweep_scale\\scale_%d_%d", i, n_P));
            
            file_name = "plots\sweep_scale\scale_sweep.gif";
            
            k_P = mean(diff(P_scale));
            
            if(k_P > 0.0)
                gamma_0 = ones(size(obj_ref.gamma))*1.0e-3;
            else
                gamma_0 = ones(size(obj_ref.gamma));
            end
            

            for idx = 1:n_P
                [obj_sca, gamma_0, res_idx, gamma_idx] = obj_ref.scaled_version(P_scale(idx), n_R_scale, normalize_freq, N_freq, aspect_set, gamma_0);
                
                gamma(:, idx) = gamma_0;
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
                [f_n(:, idx), mode_shape(:, :, idx)] = obj_sca.modal_analysis;

                subplot(2, 6, 7:9);
                cla;
                rectangle(obj_sca);
                xlim([0 6000]);
                ylim([-1 1]*1500);
                title(sprintf("Scale: %.1f kW = %.2f %% of Ref.", obj_sca.P_rated, gm_P(idx)*100.0));
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
                        yline(1.25, "color", "k");
                        ylim([1.0 2.0]);
                        yticks(1.0:0.25:2.0);
                        b(1).FaceColor = plot_prop1{end};
                        b(2).FaceColor = plot_prop2{end};
                    end
                end

                set(get(fig_axes(1), "ylabel"), "string", "$S_H$, [-]"        , latex_setting{:});
                set(fig_axes(1)  , 'xticklabel', ["s_1", "p_{1i}", "s_2", "p_{2i}", "P_3", "W_3"]);
                
                set(fig_axes, font_setting{:});
                savefig(gcf, fig_name(idx));
                print(fig_name(idx), '-dpng');
                saveasGIF(file_name, idx);
            end
            
            close all;
            
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
            
            [f_n, mode_shape] = obj_sca.modal_analysis;
            
        end
        
        function [f_n, mode_shape] = scale_resonances(obj, N, normalize, gamma_P, gamma_n, gamma, aspect)
            %SCALED_RESONANCES returns the N first resonances and mode
            % shapes of a scaled Drivetrain object. The resonances can be
            % normalized or not.
            %
            
            obj_sca = scale_aspect(obj, gamma_P, gamma_n, gamma, aspect);
            
            [f_n, mode_shape] = obj_sca.resonances(N, normalize);
            
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
