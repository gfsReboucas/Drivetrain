classdef Shaft
    %SHAFT This class implements some geometrical concepts and parameters 
    % for cylindrical shafts.
    %
    % References:
    % [1]  Budynas, R., Nisbett, J. (2015). Shigley's Mechanical
    % Engineering Design. 10th ed. New York: McGraw-Hill
    % [2] Neto M. A., Amaro A., Roseiro L., Cirne J., Leal R. (2015) Finite 
    % Element Method for Beams. In: Engineering Computation of Structures: 
    % The Finite Element Method. Springer, Cham 
    % https://doi.org/10.1007/978-3-319-17710-6_4
    %
    % written by:
    % Geraldo Rebouças
    % - geraldo.reboucas@ntnu.no OR
    % - gfs.reboucas@gmail.com
    %
    % Postdoctoral Fellow at:
    % Norwegian University of Science and Technology, NTNU
    % Department of Marine Technology, IMT
    % Marine System Dynamics and Vibration Lab, MD Lab
    % https://www.ntnu.edu/imt/lab/md-lab
    %
    
    properties
        d        (1, :) {mustBeNumeric, mustBeFinite, mustBePositive} = 1.0;   % [mm], Diameter
        L        (1, :) {mustBeNumeric, mustBeFinite, mustBePositive} = 100.0; % [mm], Length
        bearing  (1, :) Bearing = Bearing();                                   % [-],  Bearing array
        material (1, :) Material = Material();                                 % [-],  material
    end
    
    properties(Dependent)
        slender_ratio;    % [-],   Slenderness ratio
        volume;           % [m^3], Volume
        mass;             % [kg],  Mass
        A;                % [m^2], Cross section area
        I_x;              % [m^4], Area moment of inertia (rot. axis)
        I_y;              % [m^4], Area moment of inertia
        I_z;              % [m^4], Area moment of inertia
        J_x;              % [kg-m^2], Mass moment of inertia (rot. axis)
        J_y;              % [kg-m^2], Mass moment of inertia
        J_z;              % [kg-m^2], Mass moment of inertia
        r_g;              % [m], Radius of gyration
    end
    
    methods
        function obj = Shaft(dd, LL, brg, mat)
            if(nargin == 0) % LSS
                dd = 700.0;
                LL = 2.0e3;
                brg = Bearing();
                mat = Material();
            end
            
            obj.d = dd;
            obj.L = LL;
            obj.bearing = brg;
            obj.material = mat;
%             obj.F = FF;
%             obj.M = MM;
        end
        
        function tab = disp(obj)
            tab_str = {"Diameter",                           "d",    "mm",     obj.d;
                       "Length",                             "L",    "mm",     obj.L;
                       "Slenderness ratio"                   "s",    "-",      obj.slender_ratio;
                       "Mass",                               "m",    "kg",     obj.mass;
                       "Area moment of inertia (rot. axis)", "I_xx", "m^4",    obj.I_x;
                       "Mass moment of inertia (rot. axis)", "J_xx", "kg-m^2", obj.J_x;
                        };
%                        "Applied force (rot. axis)",          "F_x",  "N",      obj.F(1);
%                        "Applied force",                      "F_y",  "N",      obj.F(2);
%                        "Applied force",                      "F_z",  "N",      obj.F(3);
%                        "Applied torque (rot. axis)",         "M_x",  "N-m",    obj.M(1);
%                        "Applied torque",                     "M_y",  "N-m",    obj.M(2);
%                        "Applied torque",                     "M_z",  "N-m",    obj.M(3);
%                        "Safety factor",                      "S_F",  "-",    obj.S_H;
            
            Parameter = tab_str(:,1);
            Symbol    = tab_str(:,2);
            Unit      = tab_str(:,3);
            Value     = tab_str(:,4);

            tab = table(Parameter, Symbol, Value, Unit, ...
                        'variableNames', ["Parameter", "Symbol", "Value", "Unit"]);
            
            if(nargout == 0)
                disp(tab);
                clear tab;
            end
        end
        
        function h = rectangle(obj, varargin)
            if(nargin == 1)
                C = zeros(2, 1);
                plot_prop = {[1.0 0.0 0.0], "edgeColor", "k", "lineStyle", "-" , "faceColor", [1.0 0.0 0.0]};
            elseif(nargin > 1)
                C = varargin{1};
                plot_prop = varargin(2:end);
            else
                error("prog:input", "Too many variables.");
            end
            
            X = 0.5*obj.L*[1 -1 -1  1] + C(1);
            Y = 0.5*obj.d*[1  1 -1 -1] + C(2);
            h = fill(X, Y, plot_prop{:});

            axis equal;
            box on;
        end
    end
    
    %% Calculations:
    methods
        function k = stiffness(obj, option)
            E  = obj.material.E*1.0e6; % [Pa], Modulus of elasticity
            G  = obj.material.G*1.0e6; % [Pa], Shear modulus
            LL = obj.L*1.0e-3;
            
            switch option
                case "axial"
                    k = E*obj.A/LL;
                case "torsional"
                    k = G*obj.I_x/LL;
                case "bending"
                    k = E*obj.I_y/LL^3;
                case "full"
                    k = [obj.stiffness("axial");
                         obj.stiffness("torsional");
                         obj.stiffness("bending")];
                otherwise
                    error("prog:input", "Option [%s] is NOT valid.", option);
            end
        end
        
        function obj_sca = scaled_Shaft(obj_ref, gamma)
            
            if(numel(gamma) ~= 2) 
                error("Scaling factor gamma should have two elements.");
            end
            
            gamma_d = gamma(1, :);      gamma_L = gamma(2, :);
            
            obj_sca = Shaft(obj_ref.d*gamma_d, ...
                            obj_ref.L*gamma_L);

        end
        
        function sca = scaled_version(obj, T_scale, option)
            
            switch(option)
                case "critical_speed"
                    
                case "axial_nat_freq"
                case "torsion_nat_freq"
                case "bending_nat_freq"
                case "all_nat_freqs"
                otherwise
                    error("prog:input", "Option [%s] is NOT valid.", option);
            end
        end
        
        function f_1 = critical_speed(obj)
            % [Hz], The shaft's first critical speed
            
            E          = 206.0e9; % [N/m^2],  Young's modulus
            LL = obj.L*1.0e-3;
            
            omega_1 = sqrt(E*obj.I_y/(obj.mass/LL))*(pi/LL)^2;
            f_1 = omega_1/(2.0*pi);
        end
        
        function [n, n_y] = safety_factors(obj, K_f, K_fs, T_m)
            %SHIGLEY Calculates the safety factors for fatigue and yielding
            % for a circular shaft according to [1].
            % Input parameters:
            % - S_ut: Tensile strength, [MPa]
            % - S_y:  Yield strength, [MPa]
            % - K_f:  Fatigue stress-concentration factor for bending
            % - K_fs: Fatigue stress-concentration factor for torsion
            % - T_m: Midrange torque, [N-m]
            %
            % Assumptions:
            % - Solid shaft with round cross section.
            % - Axial loads are neglected.
            % - Using the distortion energy failure theory
            %
            % S_f:  Finite life strength 
            
            dd = obj.d;
            S_ut = obj.material.S_ut;
            S_y = obj.material.S_y;
            
            % Endurance limit
            if(S_ut <= 1400) % [MPa]
                S_e_prime = 0.5*S_ut;
            else
                S_e_prime = 700.0;
            end
            
            % Surface factor:
%             k_a =  1.58*power(S_ut, -0.085); % Ground
            k_a =  4.51*power(S_ut, -0.265); % Machined or cold-drawn
%             k_a =  57.7*power(S_ut, -0.718); % Hot-rolled
%             k_a = 272.0*power(S_ut, -0.995); % As-forged
            
            % Size factor:
            if((0.11*25.4 <= dd) || (dd <= 2.0*25.4)) % [mm]
                k_b = power(dd/7.62, -0.107);
            elseif((2.0*25.4 <= dd) || (dd <= 10.0*25.4))
                k_b = 1.51*power(dd, -0.157);
            end
            
            % Loading factor:
%             k_c = 1.0;  % Bending
%             k_c = 0.85; % Axial
            k_c = 0.59; % Torsion
            
            % Temperature factor: 70 <= T_F <= 1000 Fahrenheit
            T_F = 0.0;
            k_d = 0.975 + (0.432e-3)*T_F - (0.115e-5)*T_F^2 + (0.104e-8)*T_F^3 - (0.595e-12)*T_F^4;
            
            % Reliability factor:
%             z_a = 1.288; % 90%
            z_a = 1.645; % 95%
%             z_a = 2.326; % 99%
            k_e = 1.0 - 0.08*z_a;
            
            % Miscellaneous-Effects factor:
            k_f = 1.0;
            
            % Endurance limit at the critical location of a machine part in
            % the geometry and condition of use:
            S_e = k_a*k_b*k_c*k_d*k_e*k_f*S_e_prime;
            
            M_m = 0.0;                   % Midrange bending moment
            M_a = 0.0; % Alternating bending moment
%             T_m = max(obj.M([1 4]));     % Midrange torque
            T_a = 0.0;                   % Alternating torque
            
            dd = obj.d*1.0e-3;
            S_e = S_e*1.0e6;
            S_y = S_y*1.0e6;
            % Distortion energy theory and ASME fatigue criteria:
            val = (16.0/(pi*dd^3))*sqrt(4.0*(K_f*M_a/S_e)^2 + 3.0*(K_fs*T_a/S_e)^2 + ...
                                      4.0*(K_f*M_m/S_y)^2 + 3.0*(K_fs*T_m/S_y)^2);
            n = 1.0/val;
            % n is proportional to d^3 T^-1
            
            % Yielding factor of safety:
            sigma_Max = sqrt(3.0*(16.0*T_m/(pi*dd^3)));
            
            n_y = S_y/sigma_Max;
            
        end
        
        function M = inertia_matrix(obj, option)
            switch option
                case "axial"
                    M = eye(2)*2;
                    M(2, 1) = 1.0;      M(1, 2) = 1.0;

                    M = M*(obj.mass/6.0);

                case "torsional"
                    M = eye(2)*2;
                    M(2, 1) = 1.0;      M(1, 2) = 1.0;

                    M = M*(obj.mass*obj.I_x)/(6.0*obj.A);

                case "bending"
                    LL = obj.L*1.0e-3;
                    
                    M_d = diag([156.0 4.0*LL^2 156.0 4.0*LL^2]);
                    M_u = zeros(4);
                    M_u(1, 2) = 22.0*LL; M_u(1, 3) = 54.0;    M_u(1, 4) = -13.0*LL;
                                         M_u(2, 3) = 13.0*LL; M_u(2, 4) = - 3.0*LL^2;
                                                              M_u(3, 4) = -22.0*LL;
                    
                    M = M_u + M_u' + M_d;
                    M = M*(obj.mass/420.0);
                    
                case "full"
                    % +----------+-----------------------------------------------------------------------+
                    % |          |                             Full beam (3D)                            |
                    % +----------+-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+
                    % | DOF name | x_1 | y_1 | z_1 | a_1 | b_1 | g_1 | x_2 | y_2 | z_2 | a_2 | b_2 | g_2 |
                    % +----------+-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+
                    % |    v_j   |  1  |  2  |  3  |  4  |  5  |  6  |  7  |  8  |  9  |  10 |  11 |  12 |
                    % +----------+-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+
                    % |          |                                Initial                                |
                    % +----------+-----------+-----------+-----------------------+-----------------------+
                    % |   Load   |   Axial   | Torsional |       Bending xy      |       Bending xz      |
                    % +----------+-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+
                    % | DOF name | x_1 | x_2 | a_1 | a_2 | y_1 | g_1 | y_2 | g_2 | z_1 | b_1 | z_2 | b_2 |
                    % +----------+-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+
                    % |    u_i   |  1  |  2  |  3  |  4  |  5  |  6  |  7  |  8  |  9  |  10 |  11 |  12 |
                    % +----------+-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+
                    % |  will be |  1  |  7  |  4  |  10 |  2  |  6  |  8  |  12 |  3  |  5  |  9  |  11 |
                    % +----------+-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+
                    %
                    % Initial: M u'' + K u = f
                    % Coordinate transformation: u = R v
                    % Full beam (3D): (R^T M R) v'' + (R^T K R) v = R^T f
                    %
                    M_a = obj.inertia_matrix("axial");
                    M_t = obj.inertia_matrix("torsional");
                    M_b = obj.inertia_matrix("bending");
                    
                    M = blkdiag(M_a, ... % [1: x_1      2: x_2]
                                M_t, ... % [3: alpha_1  4: alpha_2]
                                M_b, ... % [5: y_1      6: gamma_1   7: y_2  8: gamma_2]
                                M_b);    % [9: z_1     10: beta_1   11: z_2 12: beta_2]
                    
                    R = zeros(12);
                    R( 1,  1) =  1;     R( 2,  7) = 1;     R( 3,  4) =  1;
                    R( 4, 10) =  1;     R( 5,  2) = 1;     R( 6,  6) =  1;
                    R( 7,  8) =  1;     R( 8, 12) = 1;     R( 9,  3) =  1;
                    R(10,  5) = -1;     R(11,  9) = 1;     R(12, 11) = -1;
                    
                    M = R' * M * R;

                case 'Lin_Parker_99'
                    % +----------+-----------------------------------------------------------------------+
                    % |          |                             Lin_Parker_99                             |
                    % +----------+-----------+-----------+-----------+-----------+-----------+-----------+
                    % | DOF name |    y_1    |    z_1    |    a_1    |    y_2    |    z_2    |    a_2    |
                    % +----------+-----------+-----------+-----------+-----------+-----------+-----------+
                    % |    v_j   |     1     |     2     |     3     |     4     |     5     |     6     |
                    % +----------+-----------+-----------+-----------+-----------+-----------+-----------+
                    % |          |                             Full beam (3D)                            |
                    % +----------+-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+
                    % | DOF name | x_1 | y_1 | z_1 | a_1 | b_1 | g_1 | x_2 | y_2 | z_2 | a_2 | b_2 | g_2 |
                    % +----------+-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+
                    % |    u_j   |  1  |  2  |  3  |  4  |  5  |  6  |  7  |  8  |  9  |  10 |  11 |  12 |
                    % +----------+-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+
                    % |  will be |  -  |  1  |  2  |  3  |  -  |  -  |  -  |  4  |  5  |  6  |  -  |  -  |
                    % +----------+-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+
                    %
                    % Initial: M u'' + K u = f
                    % Coordinate transformation: u = R v
                    % Lin_Parker_99: (R^T M R) v'' + (R^T K R) v = R^T f
                    %
                    
                    M = obj.inertia_matrix('full');
                    
                    R = zeros(12, 6);
                    R(2, 1) = 1;    R(3, 2) = 1;    R( 4, 3) = 1;
                    R(8, 4) = 1;    R(9, 5) = 1;    R(10, 6) = 1;
                    
                    M = R' * M * R;
                    
                otherwise
                    error("prog:input", "Option [%s] is NOT valid.", option);
            end
        end
        
        function K = stiffness_matrix(obj, option)
            E  = obj.material.E*1.0e6; % [Pa], Modulus of elasticity
            LL = obj.L*1.0e-3;
            
            switch option
                case "axial"
                    % [x_1 x_2]
                    K = eye(2);
                    K(1, 2) = -1.0;      K(2, 1) = K(1, 2);
                    
                    K = K*(E*obj.A/LL);

                case "torsional"
                    % [alpha_1 alpha_2]
                    G  = obj.material.G*1.0e6; % [Pa], Shear modulus
                    
                    K = eye(2);
                    K(1, 2) = -1.0;      K(2, 1) = K(1, 2);

                    K = K*(G*obj.I_x/LL);

                case "bending"
                    % [y_1 gamma_1 y_2 gamma_2] OR [z_1 beta_1 z_2 beta_2] 
                    LL = obj.L*1.0e-3;
                    K_d = diag([12.0 4.0*LL^2 12.0 4.0*LL^2]);
                    
                    K_u = zeros(4);
                    K_u(1, 2) = 6.0*LL; K_u(1, 3) = -12.0;    K_u(1, 4) =  6.0*LL;
                                        K_u(2, 3) =  -6.0*LL; K_u(2, 4) =  2.0*LL^2;
                                                              K_u(3, 4) = -6.0*LL;
                    K = K_u + K_u' + K_d;
                    K = K*(E*obj.I_y/LL^3);
                    
                case "full"
                    
                    K_a = obj.stiffness_matrix("axial");     % [x_1 x_2]
                    K_t = obj.stiffness_matrix("torsional"); % [alpha_1 alpha_2]
                    K_b = obj.stiffness_matrix("bending");   % [y_1 beta_1 y_2 beta_2]
                    
                    K = blkdiag(K_a, ... % [1: x_1      2: x_2]
                                K_t, ... % [3: alpha_1  4: alpha_2]
                                K_b, ... % [5: y_1      6: gamma_1   7: y_2  8: gamma_2]
                                K_b);    % [9: z_1     10: beta_1   11: z_2 12: beta_2]
                    
                    R = zeros(12);
                    R( 1,  1) =  1;     R( 2,  7) = 1;     R( 3,  4) =  1;
                    R( 4, 10) =  1;     R( 5,  2) = 1;     R( 6,  6) =  1;
                    R( 7,  8) =  1;     R( 8, 12) = 1;     R( 9,  3) =  1;
                    R(10,  5) = -1;     R(11,  9) = 1;     R(12, 11) = -1;
                    
                    K = R' * K * R;

                case 'Lin_Parker_99'
                    
                    K = obj.stiffness_matrix('full');
                    
                    R = zeros(12, 6);
                    R(2, 1) = 1;    R(3, 2) = 1;    R( 4, 3) = 1;
                    R(8, 4) = 1;    R(9, 5) = 1;    R(10, 6) = 1;
                    
                    K = R' * K * R;
                    
                otherwise
                    error("prog:input", "Option [%s] is NOT valid.", option);
            end
        end
    end
    
    methods(Static)
        function flag = test()
            %TEST
            % For a shaft with d = 50 mm and L = 1000 mm:
            %     Reference    Calculated     Rel_Diff      has_Problem
            %     _________    __________    ___________    ___________
            %
            %           0             0              NaN       "no"
            %           0             0              NaN       "no"
            %           0             0              NaN       "no"
            %           0             0              NaN       "no"
            %           0             0              NaN       "no"
            %           0             0              NaN       "no"
            %      273.81        273.81      -1.4532e-13       "no"
            %      273.81        273.81        2.076e-13       "no"
            %      935.24        935.24      -1.2156e-14       "no"
            %      935.24        935.24                0       "no"
            %      1753.8        1753.8                0       "no"
            %      2827.9        2827.9      -1.6081e-14       "no"
            %
            
            dd = 50.0;
            LL = 1.0e3;
            shaft = Shaft(dd, LL);
            
            LL = LL*1.0e-3;
            E = shaft.material.E*1.0e6;
            G = shaft.material.G*1.0e6;
            rho = shaft.material.rho*1.0e9;
            
            c_o = sqrt(E/rho);
            c_ot = sqrt(G/rho);
            
            w_na = pi*c_o/LL;
            w_nt = pi*c_ot/LL;
            
            w_nb = c_o*shaft.r_g*([4.73 7.8532]/LL).^2;
            
            fn_ana = sort([w_na; w_nt; w_nb'])/(2*pi);
            
            K_f = shaft.stiffness_matrix("full");
            M_f = shaft.inertia_matrix("full");
            
            K_LP = shaft.stiffness_matrix("Lin_Parker_99");
            M_LP = shaft.inertia_matrix("Lin_Parker_99");
            
            import matlab.mock.TestCase;
            dyn_for_TC = TestCase.forInteractiveUse;

            id = ["ISO_6336:KV", ...
                  "ISO_6336:SF", ...
                  "ISO_6336:KS", ...
                  "Dynamic_Formulation:imag", ...
                  "Dynamic_Formulation:RB"];
            
            for idx = 1:length(id)
                warning("off", id(idx));
            end
            
            [dyn_for_mock, ~] = createMock(dyn_for_TC, ?Dynamic_Formulation, 'ConstructorInputs', {NREL_5MW()});
            dyn_for_mock.K = K_f;
            dyn_for_mock.M = M_f;
            fn_ref = dyn_for_mock.modal_analysis();
            fn_ref = sort(fn_ref);
            
            name = ["axial", "torsional", "bending", "bending"];
            fn_test = [];
            for idx = 1:length(name)
                dyn_for_mock.K = shaft.stiffness_matrix(name{idx});
                dyn_for_mock.M = shaft.inertia_matrix(name{idx});
                fn_test = [fn_test;
                           dyn_for_mock.modal_analysis()];
            end
            fn_test = sort(fn_test);
            
            for idx = 1:length(id)
                warning("on", id(idx));
            end
            
            diff_fn = 100*(fn_ref - fn_test)./fn_ref;
            has_problem = repmat("no", length(fn_ref), 1);

            idx = abs(diff_fn) > 1.0;
            has_problem(idx) = "[YES]";
            
            table(fn_ref, fn_test, diff_fn, has_problem, 'variableNames', ...
                  {'Reference', 'Calculated', 'Rel_Diff', 'has_Problem'})

            flag = any(idx);
            
        end
    end
    
    %% Get methods:
    methods
        function val = get.slender_ratio(obj)
            % [-],    Slenderness ratio
            val = 2.0*obj.L/obj.d;
        end
        
        function val = get.A(obj)
            % [m^2], Cross section area
            val = pi*(obj.d/2.0)^2;
            val = val*1.0e-6;
        end
        
        function val = get.volume(obj)
            % [m^3], Volume
            val = obj.A*obj.L*1.0e-3;
        end
        
        function val = get.mass(obj)
            % [kg],     Mass
            rho = obj.material.rho*1.0e9; % [kg/m^3], Density
            val = rho*obj.volume;
        end
        
        function val = get.I_x(obj)
            % [m^4], Area moment of inertia (rot. axis)
            val = obj.I_y + obj.I_z;
        end
        
        function val = get.I_y(obj)
            % [m^4], Area moment of inertia
            val = (pi/4.0)*(obj.d/2.0)^4;
            val = val*1.0e-12;
        end
        
        function val = get.I_z(obj)
            % [m^4], Area moment of inertia
            val = obj.I_y;
        end
        
        function val = get.J_x(obj)
            val = (obj.mass/2.0)*(obj.d/2)^2;
            val = val*1.0e-6;
        end
        
        function val = get.J_y(obj)
            val = (obj.mass/12.0)*(3.0*(obj.d/2.0)^2 + obj.L^2);
            val = val*1.0e-6;
        end
        
        function val = get.J_z(obj)
            val = obj.J_y;
        end
        
        function val = get.r_g(obj)
            val = sqrt(obj.I_y/obj.A);
        end
    end
end
