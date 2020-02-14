classdef (HandleCompatible) Drivetrain
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
    % [3] IEC 61400-4:2012 Wind Turbines -- Part 4: Design Requirements for
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
    
    properties(Access = public)
        stage       (1, :) Gear_Set;                                                            % [-],      gearbox stages
        P_rated     (1, :)          {mustBeNumeric, mustBeFinite, mustBePositive} = 5.0e3;      % [kW],     Rated power
        n_rotor     (1, :)          {mustBeNumeric, mustBeFinite, mustBePositive} = 12.1;       % [1/min.], Rated rotor speed
        main_shaft  (1, :) Shaft                                                  = Shaft;      % [-],      Input Shaft
        m_Rotor     (1, :)          {mustBeNumeric, mustBeFinite, mustBePositive} = 110.0e3;    % [kg],     Rotor mass according to [3]
        J_Rotor     (1, :)          {mustBeNumeric, mustBeFinite, mustBePositive} = 57231535.0; % [kg-m^2], Rotor mass moment of inertia according to [6]
    end
    
    properties(SetAccess = private)
        m_Gen       (1, :)          {mustBeNumeric, mustBeFinite, mustBePositive} = 1900.0;     % [kg],     Generator mass according to [4]
        J_Gen       (1, :)          {mustBeNumeric, mustBeFinite, mustBePositive} = 534.116;    % [kg-m^2], Generator mass moment of inertia [4]
        N_stg       (1, 1)          {mustBeNumeric, mustBeFinite, mustBePositive} = 3;          % [-],      Number of stages
        % to store the values of some dependent variables:
        S_shaft_val (1, 4)          {mustBeNumeric, mustBeFinite, mustBePositive} = 1.0;    % [-],      Safey factor for the shafts
        S_H_val     (1, 6)          {mustBeNumeric, mustBeFinite, mustBePositive} = 1.25;   % [-],      Safety factor for surface durability (against pitting)
    end
    
    properties
        dynamic_model (1, :) string {mustBeMember(dynamic_model, ["Thomson_ToV", ...
                                                                  "Kahraman_1994", ...
                                                                  "Lin_Parker_1999"])} = "Thomson_ToV"; % which dynamic model should be used to perform modal analysis on the Drivetrain.
    end
    
    properties(Dependent)
        T_out;   % [N-m],    Output torque for each stage
        n_out;   % [1/min.], Output speed  for each stage
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
            
            N_DOF = 2;
            
            for idx = 1:obj.N_stg
                if(strcmp(obj.stage(idx).configuration, "parallel"))
                    tmp = obj.stage(idx).N_p + 1;
                elseif(strcmp(obj.stage(idx).configuration, "planetary"))
                    tmp = obj.stage(idx).N_p + 2;
                end
                
                N_DOF = N_DOF + tmp;
            end
            
            M = zeros(N_DOF, N_DOF);
            K = zeros(N_DOF, N_DOF);
            
            M(1, 1) = J_R;              M(end, end) = J_G;
            
            M_tmp = obj.main_shaft.inertia_matrix("torsional");
            K_tmp = obj.main_shaft.stiffness_matrix("torsional");
            
            M(1:2, 1:2) = M(1:2, 1:2) + M_tmp;        K(1:2, 1:2) = K(1:2, 1:2) + K_tmp;
            
            for idx = 1:obj.N_stg
                [M_tmp, K_tmp] = obj.stage(idx).Kahraman_1994;
                
                jdx = 5*idx - 3;
                
                if(strcmp(obj.stage(idx).configuration, "parallel"))
                    kdx = jdx:(jdx + 2);
                elseif(strcmp(obj.stage(idx).configuration, "planetary"))
                    kdx = jdx:(jdx + 5);
                end
                
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
