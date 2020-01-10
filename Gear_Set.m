classdef Gear_Set < Gear
    % This class implements SOME procedures for the calculation of the load
    % capacity of cylindrical involute gears with external or internal
    % teeth. Specifically the calculation of contact stresses for the
    % assessment of the surface durability of cylindrical gears.
    % In a planetary gear there are two different gear pairs:
    % (1) sun-planet;
    % (2) planet-ring;
    %
    % References:
    % [1] ISO 6336-1:2006 Calculation of load capacity of spur and helical
    % gears -- Part 1: Basic principles, introduction and general influence
    % factors
    % [2] ISO 6336-2:2006 Calculation of load capacity of spur and helical
    % gears -- Part 2: Calculation of surface durability (pitting)
    % [3] ISO/TR 6336-30:2017 Calculation of load capacity of spur and
    % helical gears -- Calculation examples for the application of ISO 6336
    % parts 1, 2, 3, 5
    % [4] Nejad, A. R., Guo, Y., Gao, Z., Moan, T. (2016). Development of a
    % 5 MW reference gearbox for offshore wind turbines. Wind Energy. 
    % https://doi.org/10.1002/we.1884
    % [5] IEC 61400-4:2012 Wind Turbines -- Part 4: Design Requirements for
    % wind turbine gearboxes
    % [6] ISO 21771:2007 Gears -- Cylindrical involute gears and gear pairs
    % -- Concepts and geometry
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
        configuration(1, :) string      {mustBeMember(configuration, ["parallel", "planetary"])} = "parallel"; % [-], Configuration of the gear set (e.g. parallel, planetary)
        N_p          (1, 1)             {mustBeInteger, mustBePositive}                          = 1;          % [-], Number of planets
        bearing      (1, :) Bearing;                                                                           % [-], Bearing array
    end
    
    properties
        shaft        (1, 1) Shaft;                                                                             % [-], Output shaft
        a_w          (1, :)             {mustBeFinite,  mustBePositive}                          = 13;         % [mm], Center distance
    end
    
    properties(Dependent)
        u;             % [-],         Gear ratio
        alpha_wt;      % [deg.],      Working transverse pressure angle
        eps_alpha;     % [-],         Transverse contact ratio
        eps_beta;      % [-],         Overlap ratio
        eps_gamma;     % [-],         Total contact ratio:
        cprime;        % [N/(mm-um)], Maximum single stiffness of a tooth pair
        cprime_th;     % [N/(mm-um)], Theoretical single stiffness
        c_gamma;       % [N/(mm-um)], Mean value of mesh stiffness per unit face witdh
        c_gamma_alpha; % [N/(mm-um)], Mean value of mesh stiffness per unit face witdh (used for K_v, K_Halpha, K_Falpha)
        c_gamma_beta;  % [N/(mm-um)], Mean value of mesh stiffness per unit face witdh (used for K_Hbeta, K_Fbeta)
        k_mesh;        % [N/m],       Mean value of mesh stiffness
        carrier;       % [-],         Planet carrier
    end
    
    methods
        function obj = Gear_Set(configuration, m_n, alpha_n, z, b, x, beta, k, bore_R, N_p, a_w, rack_type, bear, sha)
            if(nargin == 0)
                configuration = "parallel";
                m_n = 1.0;
                alpha_n = 20.0;
                z = 13*[1 1];
                b = 13.0;
                x = 0.0*[1 1];
                beta = 0.0;
                k = 0.0*[1 1];
                bore_R = 0.5*[1 1];
                a_w = m_n*z(1);
                rack_type = "A";
                bear = [Bearing, Bearing];
                sha = Shaft;
            end
            
            if(length(z) < 2)
                error("prog:input", "There should be at least two gears.");
            elseif(length(z) == 3)
                z(3) = -abs(z(3)); % because the ring is an internal gear
            end
            
            if((length(z) ~= length(x)) && (length(x) ~= length(k)) && (length(k) ~= length(bore_R)))
                error("prog:input", "The lengths of z, x, k and bore ratio should be the equal.");
            end
            
            if(std(m_n) ~= 0.0)
                error("prog:input", "Normal modules m_n should be equal for all gears.");
            elseif(std(alpha_n) ~= 0.0)
                error("prog:input", "Pressure angles alpha_n should be equal for all gears.");
            elseif(std(beta) ~= 0.0)
                error("prog:input", "Helix angles beta should be equal for all gears.");
            elseif(std(b) ~= 0.0)
                error("prog:input", "Face width b should be equal for all gears.");
            end
            
            obj@Gear(m_n, alpha_n, rack_type, z, b, x, beta, k, bore_R);
            
            if(strcmp(configuration, "planetary"))
                obj.configuration = configuration;
%                 [sun, planet, ring]  =  [1, 2, 3]
                obj.N_p    = N_p;
            elseif(strcmp(configuration, "parallel"))
                obj.configuration = configuration;
                obj.N_p    = 1;
            else
                error("prog:input", "Configuration [%s] is NOT defined.", configuration)
            end
            
            obj.a_w = a_w;
            obj.bearing = bear;
            obj.shaft = sha;
        end
        
        function tab = disp(obj)
            tmp_vec = NaN(size(obj.z));
            tmp_vec(2) = 1;
            tab_set = {"Gear Ratio",                            "u",       "-",      obj.u  *tmp_vec;
                       "Number of elements"                     "p",       "-",      obj.N_p*tmp_vec;
                       "Normal module",                         "m_n",     "mm",     obj.m_n*tmp_vec;
                       "Normal pressure angle",                 "alpha_n", "deg.",   obj.alpha_n*tmp_vec;
                       "Helix angle",                           "beta",    "deg.",   obj.beta*tmp_vec;
                       "Face width",                            "b",       "mm",     obj.b*tmp_vec;
                       "Center distance",                       "a_w",     "mm",     obj.a_w*tmp_vec;
                       "Number of teeth",                       "z",       "-",      obj.z;
                       "Profile shift coefficient",             "x",       "-",      obj.x;
                       "Reference diameter",                    "d",       "mm",     obj.d;
                       "Tip diameter",                          "d_a",     "mm",     obj.d_a;
                       "Root diameter",                         "d_f",     "mm",     obj.d_f;
                       "Mass",                                  "m",       "kg",     obj.mass;
                       "Mass moment of inertia (x axis, rot.)", "J_x",     "kg-m^2", obj.J_x;
                       "Mass moment of inertia (x axis)",       "J_y",     "kg-m^2", obj.J_y;
                       "Mass moment of inertia (x axis)",       "J_z",     "kg-m^2", obj.J_z;
                       "Bearing Names"                          "-+-",     "-",      tmp_vec;
                       "Bearing Types"                          "-+-",     "-",      tmp_vec;
                       };

            Parameter = tab_set(:,1);
            Symbol    = tab_set(:,2);
            Unit      = tab_set(:,3);
            Value     = tab_set(:,4);
            
            Value = cell2mat(Value);
            Value = num2cell(Value);
            
            Value(cellfun(@isnan,Value)) = {"-+-"}; %#ok<STRSCALR>
            
            if(strcmp(obj.configuration, "parallel"))
                for idx = 1:length(obj.z)
                    jdx = 3*idx - 2;
                    Value{end - 1, idx} = join(obj.bearing(jdx:jdx + 2).name, " / ");
                    Value{end    , idx} = join(obj.bearing(jdx:jdx + 2).type, " / ");
                end

                v_pinion = Value(:,1);
                v_wheel  = Value(:,2);
                
                tab = table(Parameter, Symbol, v_pinion, v_wheel, Unit, ...
                            'variableNames', ["Parameter", "Symbol", "Pinion", "Wheel", "Unit"]);
            elseif(strcmp(obj.configuration, "planetary"))
                Value{end - 1, 1} = obj.bearing(1).name;
                Value{end    , 1} = obj.bearing(1).type;
                
                Value{end - 1, 2} = join(obj.bearing(2:3).name, " / ");
                Value{end    , 2} = join(obj.bearing(2:3).type, " / ");
                
                Value{end - 1, 3} = obj.bearing(4).name;
                Value{end    , 3} = obj.bearing(4).type;
                
                v_sun = Value(:,1);
                v_pla = Value(:,2);
                v_rng = Value(:,3);
                v_car = {"-+-"; "-+-"; "-+-"; ... 
                         "-+-"; "-+-"; "-+-"; ...
                         "-+-"; "-+-"; "-+-"; ...
                         "-+-"; obj.carrier.d_a; obj.carrier.d_f; obj.carrier.mass; obj.carrier.J_x; obj.carrier.J_y; obj.carrier.J_z; join(obj.bearing(5:6).name, " / "); join(obj.bearing(5:6).type, " / ")};
                tab = table(Parameter, Symbol, v_sun, v_pla, v_rng, v_car, Unit, ...
                            'variableNames', ["Parameter", "Symbol", "Sun", "Planet", "Ring", "Carrier", "Unit"]);

            else
                error("prog:input", "Configuration [%s] is NOT defined.", obj.configuration);
            end
            
            if(nargout == 0)
                fprintf("Gear set:\n");
                disp(tab);
                fprintf("Bearings:\n");
                obj.bearing.disp;
                fprintf("Input shaft:\n");
                obj.shaft.disp;
                clear tab;
            end
        end
        
        function plot(obj)
            addpath("\\home.ansatt.ntnu.no\geraldod\Documents\MATLAB\Plot\linspecer");
            color = linspecer(4, 'qualitative');
            
            hold on;
            axis equal;
            box on;
            
            C_p = [obj.a_w, 0.0]';

            if(strcmp(obj.configuration, "parallel"))
                plot(obj.gear(1), C_p*0.0, "lineStyle", "-" , "lineWidth", 2.0, "color", color(1, :));
                plot(obj.gear(2), C_p*1.0, "lineStyle", "-" , "lineWidth", 2.0, "color", color(2, :));
                
                legend(["Pinion", "Wheel"], "location", "best", "fontName", "Times", "fontSize", 12.0);
                
            elseif(strcmp(obj.configuration, "planetary"))
%                 subplot(1, 2, 1)
                plot(obj.gear(1), C_p*0.0, "lineStyle", "-" , "lineWidth", 2.0, "color", color(1, :));
                plot(obj.gear(2), C_p*1.0, "lineStyle", "-" , "lineWidth", 2.0, "color", color(2, :));
                plot(obj.gear(3), C_p*0.0, "lineStyle", "-" , "lineWidth", 2.0, "color", color(3, :));
                plot(obj.carrier, C_p*0.0, "lineStyle", "-" , "lineWidth", 2.0, "color", color(4, :));
                
                RotXY = @(x)[cos(x), sin(x); sin(x), cos(x)];
                
                ang = 2.0*pi/obj.N_p;
                for idx = 2:obj.N_p
                    plot(obj.gear(2), RotXY(ang*(idx - 1))*C_p, "lineStyle", "-" , "lineWidth", 2.0, "color", color(2, :));
                end
                
                legend(["Sun", "Planet", "Ring", "Carrier"], "location", "best", "fontName", "Times", "fontSize", 12.0);
            end
            
            xlabel("y");
            ylabel("z");
            
        end
        
        function plot3(obj)
            addpath("\\home.ansatt.ntnu.no\geraldod\Documents\MATLAB\Plot\linspecer");
            color = linspecer(4, 'qualitative');
            
            C_p = [obj.a_w, 0.0]';

            if(strcmp(obj.configuration, "parallel"))
                plot3(obj.gear(1), C_p*0.0, "edgeColor", "none", "lineStyle", "-" , "faceColor", color(1, :));
                plot3(obj.gear(2), C_p*1.0, "edgeColor", "none", "lineStyle", "-" , "faceColor", color(2, :));
                
            elseif(strcmp(obj.configuration, "planetary"))
                plot3(obj.gear(1), C_p*0.0, "edgeColor", "none", "lineStyle", "-" , "faceColor", color(1, :));
                plot3(obj.gear(2), C_p*1.0, "edgeColor", "none", "lineStyle", "-" , "faceColor", color(2, :));
                plot3(obj.gear(3), C_p*0.0, "edgeColor", "none", "lineStyle", "-" , "faceColor", color(3, :));
                
                RotXY = @(x)[cos(x), sin(x); sin(x), cos(x)];
                
                ang = 2.0*pi/obj.N_p;
                for idx = 2:obj.N_p
                    plot3(obj.gear(2), RotXY(ang*(idx - 1))*C_p, "edgeColor", "none", "lineStyle", "none" , "faceColor", color(2, :));
                end
                
                legend(["Sun", "Planet", "Ring", "Carrier"], "location", "best", "fontName", "Times", "fontSize", 12.0);
            end
            
        end
        
        function rectangle(obj, varargin)
            if(nargin == 1)
                C = zeros(2, 1);
            else
                C = varargin{1};
            end
            
            addpath("\\home.ansatt.ntnu.no\geraldod\Documents\MATLAB\Plot\linspecer");
            color = linspecer(6, 'qualitative');
            
            hold on;
            if(strcmp(obj.configuration, "parallel"))
                C_w = [obj.b/2.0 0.0]' + C;
                C_p = [obj.b/2.0 obj.a_w]' + C;
                C_s = C_p + [obj.b + obj.shaft.L 0.0]'./2.0;
                
                rectangle(obj.gear(1), C_p, color(1, :), "edgeColor", "k", "lineStyle", "-" , "faceColor", color(1, :));
                rectangle(obj.gear(2), C_w, color(2, :), "edgeColor", "k", "lineStyle", "-" , "faceColor", color(2, :));
                rectangle(obj.shaft,   C_s, color(5, :), "edgeColor", "k", "lineStyle", "-" , "faceColor", color(5, :));
                
%                 legend([h_p h_w h_s], ["Pinion", "Wheel", "Shaft"], "location", "best", "fontName", "Times", "fontSize", 12.0);
                
            elseif(strcmp(obj.configuration, "planetary"))
                C_c = [obj.carrier.b/2.0, 0.0]' + C;
                C_g = C_c + [obj.b + obj.carrier.b 0]'./2.0;
                C_p = C_g + [0.0 obj.a_w]';
                C_s = C_g + [obj.b + obj.shaft.L 0.0]'./2.0;

                rectangle(obj.gear(3), C_g, color(3, :), "edgeColor", "k", "lineStyle", "-" , "faceColor", color(3, :));
                rectangle(obj.gear(1), C_g, color(1, :), "edgeColor", "k", "lineStyle", "-" , "faceColor", color(1, :));
                rectangle(obj.gear(2), C_p, color(2, :), "edgeColor", "k", "lineStyle", "-" , "faceColor", color(2, :));
                rectangle(obj.carrier, C_c, color(4, :), "edgeColor", "k", "lineStyle", "-" , "faceColor", color(4, :));
                rectangle(obj.shaft,   C_s, color(5, :), "edgeColor", "k", "lineStyle", "-" , "faceColor", color(5, :));
                
%                 legend([h_g h_p h_r h_c h_s], ["Sun", "Planet", "Ring", "Carrier", "Shaft"], "location", "best", "fontName", "Times", "fontSize", 12.0);
                
            end
            hold off;
        end
    end
    
    methods(Static)
        function obj = NREL_5MW(stage)
            %NREL_5MW returns the stages of the NREL 5 MW wind turbine
            % drivetrain according to [4], Table V. The values for the tip
            % alteration coefficients were taken from KISSSoft.
            %

            alpha_n = 20.0;        % [deg.],   Pressure angle (at reference cylinder)
%             P         = 5.0e3;     % [kW],     Rated power
%             n_c1      = 12.1;      % [1/min.], Carrier speed at stage 01 (input speed)
%             L_h       = 20*365*24; % [h],      Required life
%             K_A       = 1.25;      % [-],      Application factor
%             R_ah      = 0.8;       % [um],     Maximum arithmetic mean roughness for external gears according to [5], Sec. 7.2.7.2.
%             Q         = 6;         % [-],      Maximum accuracy grade
            rack_type = "A";       % [-],      Type of the basic rack from A to D

            switch(stage)
                case 1
                    p_1    = 3;            % [-],    Number of planet gears
                    m_n1   = 45.0;         % [mm],   Normal module
                    beta_1 = 0.0;          % [deg.], Helix angle (at reference cylinder)
                    b_1    = 491.0;        % [mm],   Face width
                    a_w1   = 863.0;        % [mm],   Center distance
                    z_s1   =  19;          % [-],    Number of teeth (sun)    [WHEEL]
                    z_p1   =  17;          % [-],    Number of teeth (planet) [PINION]
                    z_r1   =  56;          % [-],    Number of teeth (ring)
                    x_s1   =  0.617;       % [-],    Profile shift coefficient (sun)
                    x_p1   =  0.802;       % [-],    Profile shift coefficient (planet)
                    x_r1   = -0.501;       % [-],    Profile shift coefficient (ring)
                    k_s1   = -10.861/m_n1; % [-],    Tip alteration coefficient (sun)
                    k_p1   = -10.861/m_n1; % [-],    Tip alteration coefficient (planet)
                    k_r1   =   0.0;        % [-],    Tip alteration coefficient (ring)
                    
                    bore_Rs1 = 80.0/171.0;
                    bore_Rp1 = 80.0/153.0;
                    bore_Rr1 = 1.2;
                    
                    z_1 = [z_s1 z_p1 z_r1];
                    x_1 = [x_s1 x_p1 x_r1];
                    k_1 = [k_s1 k_p1 k_r1];
                    bore_R1 = [bore_Rs1 bore_Rp1 bore_Rr1];
%                                 Sun,   Planet,   Ring,   Carrier
                    bearing_1 = Bearing.NREL_5MW(1);
                    shaft_1 = Shaft.NREL_5MW(1);
                    
                    obj = Gear_Set("planetary", m_n1, alpha_n, z_1, b_1, x_1, beta_1, k_1, bore_R1, p_1, a_w1, rack_type, bearing_1, shaft_1);

                case 2
                    p_2    = 3;          % [-],    Number of planet gears
                    m_n2   = 21.0;       % [mm],   Normal module
                    beta_2 = 0.0;        % [deg.], Helix angle (at reference cylinder)
                    b_2    = 550.0;      % [mm],   Face width
                    a_w2   = 584.0;      % [mm],   Center distance
                    z_s2   =  18;        % [-],    Number of teeth (sun)    [PINION]
                    z_p2   =  36;        % [-],    Number of teeth (planet) [WHEEL]
                    z_r2   =  93;        % [-],    Number of teeth (ring)
                    x_s2   = 0.389;      % [-],    Profile shift coefficient (sun)
                    x_p2   = 0.504;      % [-],    Profile shift coefficient (planet)
                    x_r2   = 0.117;      % [-],    Profile shift coefficient (ring)
                    k_s2   = -1.75/m_n2; % [-],    Tip alteration coefficient (sun)
                    k_p2   = -1.75/m_n2; % [-],    Tip alteration coefficient (planet)
                    k_r2   =   0.0;      % [-],    Tip alteration coefficient (ring)
                    
                    bore_Rs2 = 100.0/189.0;
                    bore_Rp2 = 95.0/189.0;
                    bore_Rr2 = 1.2;

                    z_2 = [z_s2 z_p2 z_r2];
                    x_2 = [x_s2 x_p2 x_r2];
                    k_2 = [k_s2 k_p2 k_r2];
                    bore_R2 = [bore_Rs2 bore_Rp2 bore_Rr2];
%                                 Sun,   Planet,   Ring,   Carrier
                    bearing_2 = Bearing.NREL_5MW(2);
                    shaft_2 = Shaft.NREL_5MW(2);
                    
                    obj = Gear_Set("planetary", m_n2, alpha_n, z_2, b_2, x_2, beta_2, k_2, bore_R2, p_2, a_w2, rack_type, bearing_2, shaft_2);

                case 3
                    p_3    = 1;            % [-],    Number of planet gears
                    m_n3   = 14.0;         % [mm],   Normal module
                    beta_3 = 10.0;         % [deg.], Helix angle (at reference cylinder)
                    b_3    = 360.0;        % [mm],   Face width
                    a_w3   = 861.0;        % [mm],   Center distance
                    z_13   =  24;          % [-],    Number of teeth (pinion)
                    z_23   =  95;          % [-],    Number of teeth (wheel)
                    x_13   =  0.480;       % [-],    Profile shift coefficient (pinion)
                    x_23   =  0.669;       % [-],    Profile shift coefficient (wheel)
                    k_13   =  -0.938/14.0; % [-],    Tip alteration coefficient (pinion)
                    k_23   =  -0.938/14.0; % [-],    Tip alteration coefficient (wheel)
                    
                    bore_R13 = 1809.0/3086.0;
                    bore_R23 = 3385.0/9143.0;

                    z_3 = [z_13 z_23];
                    x_3 = [x_13 x_23];
                    k_3 = [k_13 k_23];
                    bore_R3 = [bore_R13 bore_R23];
%                                 Pinion, Wheel
                    bearing_3 = Bearing.NREL_5MW(3);
                    shaft_3 = Shaft.NREL_5MW(3);
                    
                    obj = Gear_Set("parallel", m_n3, alpha_n, z_3, b_3, x_3, beta_3, k_3, bore_R3, p_3, a_w3, rack_type, bearing_3, shaft_3);
                    
                otherwise
                    error("prog:input", "Option [%d] is NOT valid.", stage);
            end
        end
        
        function Z_NT = interp_ZNT(N, line)
            switch line
                case 1
                    % St, V, GGG (perl. bai.), GTS (perl.), Eh, IF (when limited pitting is permitted)
                    x = [6.0e5 1.0e7 1.0e9 1.0e10];
                    y = [1.6   1.3   1.0   0.85];
                case 2
                    % St, V, GGG (perl. bai.), GTS (perl.), Eh, IF
                    x = [1.0e5 5.0e7 1.0e9 1.0e10];
                    y = [1.6   1.0   1.0   0.85];
                case 3
                    % GG, GGG (ferr.), NT (nitr.), NV (nitr.)
                    x = [1.0e5 2.0e6 1.0e10];
                    y = [1.3   1.0   0.85];
                case 4
                    % NV (nitrocar.)
                    x = [1.0e5 2.0e6 1.0e10];
                    y = [1.1   1.0   0.85];
                otherwise
                    error("prog:input", "Invalid input [%d].\nValid options are 1 to 4.\n", line);
            end

            if(N <= x(1))
                Z_NT = y(1);
            elseif(N > x(end))
                Z_NT = y(end);
            else
                Z_NT = interp1(log(x), y, log(N));
            end
        end
    end
    
    %% Calculation:
    methods
        %% Scaling:
        function [scaled_obj, gamma, RN] = scaled_version(obj_ref, gamma_0, SHref, P, n)
            fun_min = @(g)(SHref - obj_ref.pitting_factors(P, n, g));
            fun = @(p)(norm(fun_min(p))^2);
            
            const = @(g)deal(1.25 - obj_ref.pitting_factors(P, n, g), []);
            
            gamma_min = ones(1, 2)*1.0e-6;
            gamma_Max = ones(1, 2);
%             gamma_min(1) = 1.0/obj_ref.m_n;
%             gamma_Max(1) = 50.0/obj_ref.m_n;
            
            opt_solver = optimoptions("fmincon", "Display", "notify");

            id = "prog:input";
            warning("off", id);
            [gamma, RN, ~]  = fmincon(fun, gamma_0*ones(1, 2), [], [], [], [], gamma_min, gamma_Max, const, opt_solver);
            warning("on", id);
            
            scaled_obj = obj_ref.scaled_Gear_Set(gamma);
        end
        
        function scaled_obj = scaled_Gear_Set(obj, gamma)
            if(length(gamma) == 1)
                m_n_sca = Rack.module(obj.m_n*gamma, "calc", "nearest");
                b_sca   =             obj.b  *gamma;
            elseif(length(gamma) == 2)
                m_n_sca = Rack.module(obj.m_n*gamma(1), "calc", "nearest");
                b_sca   =             obj.b  *gamma(2);
            else
                error("prog:input", "Wrong size for gamma. %d ~= 1 or 3.", length(gamma));
            end
            
            a_w_sca = obj.a_w*(m_n_sca/obj.m_n);
            
            scaled_obj = obj.modify_Gear_Set(a_w_sca, b_sca, m_n_sca);
        end
        
        function mod_obj = modify_Gear_Set(obj, aw, bb, mn)
            mod_obj = Gear_Set(obj.configuration, ... % planetary or parallel
                               mn,                ... % normal module                  
                               obj.alpha_n,       ... % pressure angle
                               abs(obj.z),        ... % number of teeth
                               bb,                ... % face width
                               obj.x,             ... % profile shift coefficient
                               obj.beta,          ... % helix angle
                               obj.k,             ... % tip alteration coefficient
                               obj.bore_ratio,    ... % ratio btw bore and reference diameters
                               obj.N_p,           ... % number of planets
                               aw,                ... % center distance
                               obj.type,          ... % rack type
                               obj.bearing,       ... % bearing array
                               obj.shaft);        ... % shaft
            
            if(abs(obj.alpha_wt - mod_obj.alpha_wt) > 1.0e-3)
                mod_obj.a_w = mod_obj.find_center_distance(obj.alpha_wt);
            end
        end

        %% Dynamics:
        function [M, G, K_Omega] = dynamic_matrix(obj, option)
            switch option
                case "y-z-alpha" % "Lin_Parker_1999"
                    [M, G, K_Omega] = obj.Lin_Parker_1999;
                case "alpha" % "Kahraman_1994"
                    M       = obj.Kahraman_1994;
                    G       = zeros(size(M));
                    K_Omega = zeros(size(M));
                case "6DOF" % "Eritenel"
                    error("prog:input", " Option [%s] is NOT implemented yet.", upper(option));
                otherwise
                    error("prog:input", " Option [%s] is NOT valid.", upper(option));
            end
        end
        
        function [M, K] = Kahraman_1994(obj)
            %KAHRAMAN_1994 Returns the inertia and stiffness matrices of
            % the gear stage according to:
            % A. Kahraman, "Natural Modes of Planetary Gear Trains",
            % Journal of Sound and Vibration, vol. 173, no. 1, pp. 125-130,
            % 1994. https://doi.org/10.1006/jsvi.1994.1222.

            
            M = NaN;
            K = NaN;

            if(strcmp(obj.configuration, "parallel"))
                n = 3;
                M = zeros(n, n);                K = zeros(n, n);
                
                m_p = obj.mass(1);
                m_w = obj.mass(2);
                
                r_p = obj.d(1)*1.0e-3/2.0;
                r_w = obj.d(2)*1.0e-3/2.0;
                
                k_pw = obj.k_mesh;
                
                idx = 1:2;
                M(idx, idx) = diag([m_w*r_w^2 m_p*r_p^2]);
                
                K(idx, idx) = [r_w ^ 2 * k_pw             r_p     * k_pw * r_w;
                               r_w     * k_pw * r_p   3 * r_p ^ 2 * k_pw;];
                
            elseif(strcmp(obj.configuration, "planetary"))
                n = 6;
                M = zeros(n, n);                K = zeros(n, n);
                
                m_1 = obj.mass(1);
                m_2 = obj.mass(2);
                m_c = obj.carrier.mass;
                
                aw = obj.a_w*1.0e-3/2.0;
                r_1 = obj.d(1)*1.0e-3/2.0;
                r_2 = obj.d(2)*1.0e-3/2.0;

                % Mesh stiffness:
                sun_pla = Gear_Set("parallel",          obj.m_n, ...
                                   obj.alpha_n,         obj.z(1:2), ...
                                   obj.b,               obj.x(1:2), ...
                                   obj.beta,            obj.k(1:2), ...
                                   obj.bore_ratio(1:2), obj.N_p, ...
                                   obj.a_w,             obj.type, ...
                                   obj.bearing,         obj.shaft);

                pla_rng = Gear_Set("parallel",          obj.m_n, ...
                                   obj.alpha_n,         obj.z(2:3), ...
                                   obj.b,               obj.x(2:3), ...
                                   obj.beta,            obj.k(2:3), ...
                                   obj.bore_ratio(2:3), obj.N_p, ...
                                   obj.a_w,             obj.type, ...
                                   obj.bearing,         obj.shaft);

                k_sp = sun_pla.k_mesh;
                k_rp = pla_rng.k_mesh;
                
                idx = 1:5;
                M(idx, idx) = diag([m_c*aw^2 ones(1, 3)*m_2*r_2^2 m_1*r_1^2]);
                K(idx, idx) = [3 * aw ^ 2 * (k_rp + k_sp) aw * (k_rp - k_sp) * r_2 aw * (k_rp - k_sp) * r_2 aw * (k_rp - k_sp) * r_2 -3 * aw * k_sp * r_1;
                                   aw     * (k_rp - k_sp) * r_2 r_2 ^ 2 * (k_rp + k_sp) 0 0 r_2 * k_sp * r_1;
                                   aw     * (k_rp - k_sp) * r_2 0 r_2 ^ 2 * (k_rp + k_sp) 0 r_2 * k_sp * r_1;
                                   aw     * (k_rp - k_sp) * r_2 0 0 r_2 ^ 2 * (k_rp + k_sp) r_2 * k_sp * r_1;
                              -3 * aw * k_sp * r_1 r_2 * k_sp * r_1 r_2 * k_sp * r_1 r_2 * k_sp * r_1 r_1 ^ 2 * (3 * k_sp)];
            end
            
            idx = (n - 1):n;
            
            M(idx, idx) = M(idx, idx) + obj.shaft.inertia_matrix("torsional");
            K(idx, idx) = K(idx, idx) + obj.shaft.stiffness_matrix("torsional");
        end
        
        function [M, K, K_b, K_m, K_Omega, G] = Lin_Parker_1999(obj)
            %LIN_PARKER_1999 Returns the inertia and stiffness matrices of
            % the gear set according to:
            % J. Lin and R. Parker, "Analytical Characterization of the
            % Unique Properties of Planetary Gear Free Vibration", Journal
            % of Vibration and Acoustics, vol. 121, no. 3, pp. 316-321,
            % 1999. https://doi.org/10.1115/1.2893982
            %
            
            if(strcmp(obj.configuration, "parallel"))
                n = 9;
                M       = zeros(n, n);
                K_b     = zeros(n, n);
                K_m     = zeros(n, n);
                K_Omega = zeros(n, n);
                G       = zeros(n, n);
                
                J_p = obj.J_x(1);               J_w = obj.J_x(2);
                m_p = obj.mass(1);              m_w = obj.mass(2);
                r_p = obj.d(1)*1.0e-3/2.0;      r_w = obj.d(2)*1.0e-3/2.0;
                
                alpha_nn = obj.alpha_n;
                
                bear = obj.bearing;
                b_p = parallel_association(bear(4:6));
                b_w = parallel_association(bear(1:3));

                % [x, y, theta] = [y, z, alpha]
                k_px = b_p.K_y;     k_py = b_p.K_z;     k_pu = b_p.K_alpha;
                k_wx = b_w.K_y;     k_wy = b_w.K_z;     k_wu = b_w.K_alpha;
                
                k_p = obj.k_mesh;
                
                idx = 1:(n - 3);
                
                % Mass matrix
                M(idx, idx)   = diag([m_w m_w J_w ...
                                      m_p m_p J_p]);

                % Bearing stiffness:
                K_b(idx, idx) = diag([zeros(1,3), k_px, k_py, k_pu*r_p^2]);
                
                K_m(idx, idx) = [-k_p*cosd(alpha_nn)^2 + k_p + k_wx,  k_p*cosd(alpha_nn)*sind(alpha_nn), -r_w*k_p*sind(alpha_nn),              -k_p*sind(alpha_nn)^2, -k_p*cosd(alpha_nn)*sind(alpha_nn), -r_p*k_p*sind(alpha_nn);
                                  k_p*cosd(alpha_nn)*sind(alpha_nn),        k_wy + k_p*cosd(alpha_nn)^2, -r_w*k_p*cosd(alpha_nn), -k_p*cosd(alpha_nn)*sind(alpha_nn),              -k_p*cosd(alpha_nn)^2, -r_p*k_p*cosd(alpha_nn);
                                            -r_w*k_p*sind(alpha_nn),            -r_w*k_p*cosd(alpha_nn),      r_w^2*(k_wu + k_p),             r_w*k_p*sind(alpha_nn),             r_w*k_p*cosd(alpha_nn),             r_w*k_p*r_p;
                                              -k_p*sind(alpha_nn)^2, -k_p*cosd(alpha_nn)*sind(alpha_nn),  r_w*k_p*sind(alpha_nn),               k_p*sind(alpha_nn)^2,  k_p*cosd(alpha_nn)*sind(alpha_nn),  r_p*k_p*sind(alpha_nn);
                                 -k_p*cosd(alpha_nn)*sind(alpha_nn),              -k_p*cosd(alpha_nn)^2,  r_w*k_p*cosd(alpha_nn),  k_p*cosd(alpha_nn)*sind(alpha_nn),               k_p*cosd(alpha_nn)^2,  r_p*k_p*cosd(alpha_nn);
                                            -r_p*k_p*sind(alpha_nn),            -r_p*k_p*cosd(alpha_nn),             r_w*k_p*r_p,             r_p*k_p*sind(alpha_nn),             r_p*k_p*cosd(alpha_nn),              r_p^2*k_p];
                K_Omega(idx, idx) = diag([m_w m_w 0.0 ...
                                          m_p m_p 0.0]);
                
                % Gyroscopic matrix:
                G(idx, idx) = [0.0     -2.0*m_w  0.0      0.0      0.0      0.0;
                               2.0*m_w  0.0      0.0      0.0      0.0      0.0;
                               0.0      0.0      0.0      0.0      0.0      0.0;
                               0.0      0.0      0.0      0.0     -2.0*m_p  0.0;
                               0.0      0.0      0.0      2.0*m_p  0.0      0.0;
                               0.0      0.0      0.0      0.0      0.0      0.0];

             elseif(strcmp(obj.configuration, "planetary"))
                n = 18;
                M       = zeros(n, n);
                K_b     = zeros(n, n);
                K_m     = zeros(n, n);
                K_Omega = zeros(n, n);
                G       = zeros(n, n);

                J_s = obj.J_x(1);
                J_p = obj.J_x(2);
                J_c = obj.carrier.J_x;
                
                m_s = obj.mass(1);
                m_p = obj.mass(2);
                m_c = obj.carrier.mass;
                
                r_s = obj.d(1)*1.0e-3/2.0;
                r_p = obj.d(2)*1.0e-3/2.0;
                r_c = obj.a_w*1.0e-3;
                
                alpha_nn = obj.alpha_n;
                
                psi = 360.0*((1:3) - 1)/obj.N_p;
                
                bear = obj.bearing;
                b_s = parallel_association(bear(1:2));
                b_p = parallel_association(bear(3:4));
                b_c = parallel_association(bear(7:8));
                
                % [x, y, theta] = [y, z, alpha]
                k_sx = b_s.K_y;     k_sy = b_s.K_z;     k_su = b_s.K_alpha;
                k_px = b_p.K_y;     k_py = b_p.K_z;     k_pu = b_p.K_alpha;
                k_cx = b_c.K_y;     k_cy = b_c.K_z;     k_cu = b_c.K_alpha;

                % Mesh stiffness:
                sun_pla = Gear_Set("parallel", obj.m_n, obj.alpha_n,     obj.z(1:2) , obj.b, obj.x(1:2), obj.beta, obj.k(1:2), obj.bore_ratio(1:2), 0, obj.a_w, obj.type, obj.bearing, obj.shaft);
                pla_rng = Gear_Set("parallel", obj.m_n, obj.alpha_n, abs(obj.z(2:3)), obj.b, obj.x(2:3), obj.beta, obj.k(2:3), obj.bore_ratio(2:3), 0, obj.a_w, obj.type, obj.bearing, obj.shaft);
                
                k_s = sun_pla.k_mesh;
                k_r = pla_rng.k_mesh;
                
                idx = 1:(n - 3);
                
                % Mass matrix
                M(idx, idx) = diag([m_c m_c J_c ...
                                   [m_p m_p J_p]*repmat(eye(3), 1, obj.N_p) ...
                                    m_s m_s J_s]);

                % Bearing stiffness:
                K_b(idx, idx) = diag([k_cx k_cy k_cu*r_c^2 ...
                                      zeros(1, 9) ...
                                      k_sx k_sy k_su*r_s^2]);

                % Mesh stiffnes matrix:
                K_m(idx, idx) = [                                              3.0*k_px,                                                   0.0, -r_c*k_pu*(sind(psi(3)) + sind(psi(2)) + sind(psi(1))),                                                     -cosd(psi(1))*k_px,                                                      sind(psi(1))*k_py,                                            0.0,                                                     -cosd(psi(2))*k_px,                                                      sind(psi(2))*k_py,                                            0.0,                                                     -cosd(psi(3))*k_px,                                                      sind(psi(3))*k_py,                                            0.0,                                                                                                                                                             0.0,                                                                                                                                                             0.0,                                                                                        0;
                                                                                    0.0,                                              3.0*k_py,  r_c*k_pu*(cosd(psi(3)) + cosd(psi(2)) + cosd(psi(1))),                                                     -sind(psi(1))*k_px,                                                     -cosd(psi(1))*k_py,                                            0.0,                                                     -sind(psi(2))*k_px,                                                     -cosd(psi(2))*k_py,                                            0.0,                                                     -sind(psi(3))*k_px,                                                     -cosd(psi(3))*k_py,                                            0.0,                                                                                                                                                             0.0,                                                                                                                                                             0.0,                                                                                        0;
                                 -r_c*k_px*(sind(psi(3)) + sind(psi(2)) + sind(psi(1))), r_c*k_py*(cosd(psi(3)) + cosd(psi(2)) + cosd(psi(1))),                                         3.0*r_c^2*k_pu,                                                                    0.0,                                                              -r_c*k_py,                                            0.0,                                                                    0.0,                                                              -r_c*k_py,                                            0.0,                                                                    0.0,                                                              -r_c*k_py,                                            0.0,                                                                                                                                                             0.0,                                                                                                                                                             0.0,                                                                                        0;
                                                                     -cosd(psi(1))*k_px,                                    -sind(psi(1))*k_px,                                                    0.0,        -k_r*cosd(alpha_nn)^2 - k_s*cosd(alpha_nn)^2 + k_r + k_s + k_px, -k_r*cosd(alpha_nn)*sind(alpha_nn) + k_s*cosd(alpha_nn)*sind(alpha_nn), r_p*(-k_r*sind(alpha_nn) - k_s*sind(alpha_nn)),                                                                    0.0,                                                                    0.0,                                            0.0,                                                                    0.0,                                                                    0.0,                                            0.0,                                                                                                                    -k_s*sind(-psi(1) + alpha_nn)*sind(alpha_nn),                                                                                                                    -k_s*cosd(-psi(1) + alpha_nn)*sind(alpha_nn),                                                                  -r_s*k_s*sind(alpha_nn);
                                                                      sind(psi(1))*k_py,                                    -cosd(psi(1))*k_py,                                              -r_c*k_py, -k_r*cosd(alpha_nn)*sind(alpha_nn) + k_s*cosd(alpha_nn)*sind(alpha_nn),                     k_py + k_r*cosd(alpha_nn)^2 + k_s*cosd(alpha_nn)^2,  r_p*(k_r*cosd(alpha_nn) - k_s*cosd(alpha_nn)),                                                                    0.0,                                                                    0.0,                                            0.0,                                                                    0.0,                                                                    0.0,                                            0.0,                                                                                                                    -k_s*sind(-psi(1) + alpha_nn)*cosd(alpha_nn),                                                                                                                    -k_s*cosd(-psi(1) + alpha_nn)*cosd(alpha_nn),                                                                  -r_s*k_s*cosd(alpha_nn);
                                                                                    0.0,                                                   0.0,                                                    0.0,                         r_p*(-k_r*sind(alpha_nn) - k_s*sind(alpha_nn)),                          r_p*(k_r*cosd(alpha_nn) - k_s*cosd(alpha_nn)),                       r_p^2*(k_pu + k_r + k_s),                                                                    0.0,                                                                    0.0,                                            0.0,                                                                    0.0,                                                                    0.0,                                            0.0,                                                                                                                                r_p*k_s*sind(-psi(1) + alpha_nn),                                                                                                                                r_p*k_s*cosd(-psi(1) + alpha_nn),                                                                              r_p*k_s*r_s;
                                                                     -cosd(psi(2))*k_px,                                    -sind(psi(2))*k_px,                                                    0.0,                                                                    0.0,                                                                    0.0,                                            0.0,        -k_r*cosd(alpha_nn)^2 - k_s*cosd(alpha_nn)^2 + k_r + k_s + k_px, -k_r*cosd(alpha_nn)*sind(alpha_nn) + k_s*cosd(alpha_nn)*sind(alpha_nn), r_p*(-k_r*sind(alpha_nn) - k_s*sind(alpha_nn)),                                                                    0.0,                                                                    0.0,                                            0.0,                                                                                                                    -k_s*sind(-psi(2) + alpha_nn)*sind(alpha_nn),                                                                                                                    -k_s*cosd(-psi(2) + alpha_nn)*sind(alpha_nn),                                                                  -r_s*k_s*sind(alpha_nn);
                                                                      sind(psi(2))*k_py,                                    -cosd(psi(2))*k_py,                                              -r_c*k_py,                                                                    0.0,                                                                    0.0,                                            0.0, -k_r*cosd(alpha_nn)*sind(alpha_nn) + k_s*cosd(alpha_nn)*sind(alpha_nn),                     k_py + k_r*cosd(alpha_nn)^2 + k_s*cosd(alpha_nn)^2,  r_p*(k_r*cosd(alpha_nn) - k_s*cosd(alpha_nn)),                                                                    0.0,                                                                    0.0,                                            0.0,                                                                                                                    -k_s*sind(-psi(2) + alpha_nn)*cosd(alpha_nn),                                                                                                                    -k_s*cosd(-psi(2) + alpha_nn)*cosd(alpha_nn),                                                                  -r_s*k_s*cosd(alpha_nn);
                                                                                    0.0,                                                   0.0,                                                    0.0,                                                                    0.0,                                                                    0.0,                                            0.0,                         r_p*(-k_r*sind(alpha_nn) - k_s*sind(alpha_nn)),                          r_p*(k_r*cosd(alpha_nn) - k_s*cosd(alpha_nn)),                       r_p^2*(k_pu + k_r + k_s),                                                                    0.0,                                                                    0.0,                                            0.0,                                                                                                                                r_p*k_s*sind(-psi(2) + alpha_nn),                                                                                                                                r_p*k_s*cosd(-psi(2) + alpha_nn),                                                                              r_p*k_s*r_s;
                                                                     -cosd(psi(3))*k_px,                                    -sind(psi(3))*k_px,                                                    0.0,                                                                    0.0,                                                                    0.0,                                            0.0,                                                                    0.0,                                                                    0.0,                                            0.0,        -k_r*cosd(alpha_nn)^2 - k_s*cosd(alpha_nn)^2 + k_r + k_s + k_px, -k_r*cosd(alpha_nn)*sind(alpha_nn) + k_s*cosd(alpha_nn)*sind(alpha_nn), r_p*(-k_r*sind(alpha_nn) - k_s*sind(alpha_nn)),                                                                                                                    -k_s*sind(-psi(3) + alpha_nn)*sind(alpha_nn),                                                                                                                    -k_s*cosd(-psi(3) + alpha_nn)*sind(alpha_nn),                                                                  -r_s*k_s*sind(alpha_nn);
                                                                      sind(psi(3))*k_py,                                    -cosd(psi(3))*k_py,                                              -r_c*k_py,                                                                    0.0,                                                                    0.0,                                            0.0,                                                                    0.0,                                                                    0.0,                                            0.0, -k_r*cosd(alpha_nn)*sind(alpha_nn) + k_s*cosd(alpha_nn)*sind(alpha_nn),                     k_py + k_r*cosd(alpha_nn)^2 + k_s*cosd(alpha_nn)^2,  r_p*(k_r*cosd(alpha_nn) - k_s*cosd(alpha_nn)),                                                                                                                    -k_s*sind(-psi(3) + alpha_nn)*cosd(alpha_nn),                                                                                                                    -k_s*cosd(-psi(3) + alpha_nn)*cosd(alpha_nn),                                                                  -r_s*k_s*cosd(alpha_nn);
                                                                                    0.0,                                                   0.0,                                                    0.0,                                                                    0.0,                                                                    0.0,                                            0.0,                                                                    0.0,                                                                    0.0,                                            0.0,                         r_p*(-k_r*sind(alpha_nn) - k_s*sind(alpha_nn)),                          r_p*(k_r*cosd(alpha_nn) - k_s*cosd(alpha_nn)),                       r_p^2*(k_pu + k_r + k_s),                                                                                                                                r_p*k_s*sind(-psi(3) + alpha_nn),                                                                                                                                r_p*k_s*cosd(-psi(3) + alpha_nn),                                                                              r_p*k_s*r_s;
                                                                                    0.0,                                                   0.0,                                                    0.0,                           -k_s*sind(-psi(1) + alpha_nn)*sind(alpha_nn),                           -k_s*sind(-psi(1) + alpha_nn)*cosd(alpha_nn),               r_p*k_s*sind(-psi(1) + alpha_nn),                           -k_s*sind(-psi(2) + alpha_nn)*sind(alpha_nn),                           -k_s*sind(-psi(2) + alpha_nn)*cosd(alpha_nn),               r_p*k_s*sind(-psi(2) + alpha_nn),                           -k_s*sind(-psi(3) + alpha_nn)*sind(alpha_nn),                           -k_s*sind(-psi(3) + alpha_nn)*cosd(alpha_nn),               r_p*k_s*sind(-psi(3) + alpha_nn),                                                              k_s*(  3.0 - cosd(-psi(1) + alpha_nn)^2 - cosd(-psi(2) + alpha_nn)^2 - cosd(-psi(3) + alpha_nn)^2), k_s*(sind(-psi(1) + alpha_nn)*cosd(-psi(1) + alpha_nn) + sind(-psi(2) + alpha_nn)*cosd(-psi(2) + alpha_nn) + sind(-psi(3) + alpha_nn)*cosd(-psi(3) + alpha_nn)), r_s*k_s*(sind(-psi(1) + alpha_nn) + sind(-psi(2) + alpha_nn) + sind(-psi(3) + alpha_nn));
                                                                                    0.0,                                                   0.0,                                                    0.0,                           -k_s*cosd(-psi(1) + alpha_nn)*sind(alpha_nn),                           -k_s*cosd(-psi(1) + alpha_nn)*cosd(alpha_nn),               r_p*k_s*cosd(-psi(1) + alpha_nn),                           -k_s*cosd(-psi(2) + alpha_nn)*sind(alpha_nn),                           -k_s*cosd(-psi(2) + alpha_nn)*cosd(alpha_nn),               r_p*k_s*cosd(-psi(2) + alpha_nn),                           -k_s*cosd(-psi(3) + alpha_nn)*sind(alpha_nn),                           -k_s*cosd(-psi(3) + alpha_nn)*cosd(alpha_nn),               r_p*k_s*cosd(-psi(3) + alpha_nn), k_s*(sind(-psi(1) + alpha_nn)*cosd(-psi(1) + alpha_nn) + sind(-psi(2) + alpha_nn)*cosd(-psi(2) + alpha_nn) + sind(-psi(3) + alpha_nn)*cosd(-psi(3) + alpha_nn)),                                                                      k_s*(cosd(-psi(1) + alpha_nn)^2 + cosd(-psi(2) + alpha_nn)^2 + cosd(-psi(3) + alpha_nn)^2), r_s*k_s*(cosd(-psi(1) + alpha_nn) + cosd(-psi(2) + alpha_nn) + cosd(-psi(3) + alpha_nn));
                                                                                    0.0,                                                   0.0,                                                    0.0,                                                -r_s*k_s*sind(alpha_nn),                                                -r_s*k_s*cosd(alpha_nn),                                    r_p*k_s*r_s,                                                -r_s*k_s*sind(alpha_nn),                                                -r_s*k_s*cosd(alpha_nn),                                    r_p*k_s*r_s,                                                -r_s*k_s*sind(alpha_nn),                                                -r_s*k_s*cosd(alpha_nn),                                    r_p*k_s*r_s,                                                                        r_s*k_s*(sind(-psi(1) + alpha_nn) + sind(-psi(2) + alpha_nn) + sind(-psi(3) + alpha_nn)),                                                                        r_s*k_s*(cosd(-psi(1) + alpha_nn) + cosd(-psi(2) + alpha_nn) + cosd(-psi(3) + alpha_nn)),                                                                           3.0*r_s^2*k_s];

                % Centripetal stiffness matrix:
                K_Omega(idx, idx) = diag([m_c m_c 0.0 ...
                                         [m_p m_p 0.0]*repmat(eye(3), 1, obj.N_p) ...
                                          m_s m_s 0.0]);
                                  
                % Gyroscopic matrix:
                G(idx, idx) = [0.0     -2.0*m_c  0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0;
                               2.0*m_c  0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0;
                               0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0;
                               0.0      0.0      0.0      0.0     -2.0*m_p  0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0;
                               0.0      0.0      0.0      2.0*m_p  0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0;
                               0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0;
                               0.0      0.0      0.0      0.0      0.0      0.0      0.0     -2.0*m_p  0.0      0.0      0.0      0.0      0.0      0.0      0.0;
                               0.0      0.0      0.0      0.0      0.0      0.0      2.0*m_p  0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0;
                               0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0;
                               0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0     -2.0*m_p  0.0      0.0      0.0      0.0;
                               0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0      2.0*m_p  0.0      0.0      0.0      0.0      0.0;
                               0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0;
                               0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0     -2.0*m_s  0.0;
                               0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0      2.0*m_s  0.0      0.0;
                               0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0];
                
            end
            
            M_tmp = obj.shaft.inertia_matrix("full");
            K_tmp = obj.shaft.stiffness_matrix("full");
            
            % [ 1   2   3       4      5       6   7   8   9    10      11      12]
            % [x_1 y_1 z_1 alpha_1 beta_1 gamma_1 x_2 y_2 z_2 alpha_2 beta_2 gamma_2]
%             idx = [3 5 6 9 11 12];
            idx = [1 5 6 7 11 12];
            M_tmp(idx, :) = [];            M_tmp(:, idx) = [];
            K_tmp(idx, :) = [];            K_tmp(:, idx) = [];
            
%             M_tmp = diag(diag(M_tmp));           K_tmp = diag(diag(K_tmp));
            
            idx = (n - 5):n;
            
            M(  idx, idx) = M(  idx, idx) + M_tmp;
            K_b(idx, idx) = K_b(idx, idx) + K_tmp;
            
            K = @(Om)(K_b + K_m - K_Omega*Om^2);
            % remove the shaft participation ?
%             K_b(idx, idx) = K_b(idx, idx) - K_tmp;
        end
        
        function [M, K] = Eritenel_2011(obj)
            %ERITENEL_2011 Returns the inertia and stiffness matrices of
            % the drivetrain according to:
            % [1] T. Eritenel, "Three-Dimensional Nonlinear Dynamics and
            % Vibration Reduction of Gear Pairs and Planetary Gears",
            % Ph.D., Ohio State University, Mechanical Engineering, 2011.
            % http://rave.ohiolink.edu/etdc/view?acc_num=osu1298651902
            % [2] T. Eritenel and R. Parker, "Modal properties of
            % three-dimensional helical planetary gears", Journal of Sound
            % and Vibration, vol. 325, no. 1-2, pp. 397-420, 2009.
            % https://doi.org/10.1016/j.jsv.2009.03.002
            % [3] T. Eritenel and R. Parker, "Three-dimensional nonlinear
            % vibration of gear pairs", Journal of Sound and Vibration, 
            % vol. 331, no. 15, pp. 3628-3648, 2012.
            % https://doi.org/10.1016/j.jsv.2012.03.019
            % [4] T. Eritenel and R. Parker, "An investigation of tooth
            % mesh nonlinearity and partial contact loss in gear pairs
            % using a lumped-parameter model", Mechanism and Machine
            % Theory, vol. 56, pp. 28-51, 2012.
            % https://doi.org/10.1016/j.mechmachtheory.2012.05.002
            %
            
        end
        %% Pitting:
        function [S_H, sigma_H, K_Halpha, K_v, Z_B, Z_D, Z_H, Z_NT1, Z_NT2, Z_v, Z_eps] = Pitting_ISO(obj, P_inp, n_1, S_Hmin, L_h, Q, R_ah, K_A)
            %PITTING_ISO calculates the safety factor against pitting S_H
            % and the contact stress sigma_H on a gear set. The
            % calculations are based on the ISO 6336-2:2006.
            % Z_NT: dont bother
            % Z_v: gamma_n = 1/gamma_m_n [?]
            % Z_B/D:    function of a_w, b, m_n
            % K_v:      function of n_c, P/(b*n_c*m_n)
            % K_Halpha: function of b, m_n, P, n_c
            % Z_eps:    function of m_n, b
            % Z_H:      function of m_n, a_w

            %% Constants:
            E          = 206.0e3; % [N/mm^2],  Young's modulus
            nu         = 0.3;     % [-],       Poisson's ratio
            sigma_Hlim = 1500.0;  % [N/mm^2],  Allowable contact stress number
            rho        = 7.83e-6; % [kg/mm^3], Density
            % ISO viscosity grade: VG 220
            nu_40      = 220.0;   % [mm/s^2],  Nominal kinematic viscosity

%             Q    = 6;        % [-],  ISO accuracy grade
%             R_ah = 0.8;      % [um], Maximum arithmetic mean roughness for external gears according to [7], Sec. 7.2.7.2.
            R_zh = 6.0*R_ah; % [um], Mean peak-to-valley roughness

            %% Preparatory calculations:
            T_1 = (P_inp*1.0e3)/(n_1*pi/30.0);
            
            T_1 = abs(T_1);
            n_1 = abs(n_1);
            
            u_tmp = abs(obj.z(2)/obj.z(1));
            
            % Single pitch tolerance, according to ISO 1328-1:1995:
            f_pt = single_pitch_tol(obj, Q);
            
            % consider only (1:2) to discard the ring gear
            f_pb = max(f_pt(1:2).*cosd(obj.alpha_t)); % [um], Transverse base pitch deviation
            
            if(f_pb >= 40.0) %[um]
                y_alpha = 3.0; % [um]
            else
                y_alpha = f_pb*75.0e-3;
            end
            
            f_falp_tmp = profile_form_tol(obj, Q);

            % consider only (1:2) to discard the ring gear
            f_falpha = max(f_falp_tmp(1:2));

            % Estimated running allowance (pitch deviation):
            y_p = y_alpha;

            % Estimated running allowance (flank deviation):
            y_f = f_falpha*75.0e-3;

            f_pbeff = f_pb - y_p;

            f_falphaeff = f_falpha - y_f; % 0.925*y_f
            
            v = obj.pitch_line_vel(n_1);

            %% Pitting calculation:
            % Mesh load factor according to [7]:
            switch obj.N_p
                case 3
                    K_gamma = 1.1;
                case 4
                    K_gamma = 1.25;
                case 5
                    K_gamma = 1.35;
                case 6
                    K_gamma = 1.44;
                case 7
                    K_gamma = 1.47;
                otherwise
                    K_gamma = 1.0;
            end
            
            % Contact ratio factor:
            Z_eps = contact_ratio_factor(obj.beta, obj.eps_alpha, obj.eps_beta);

            % [N], Nominal tangential load
            F_t = 2.0e3*(T_1/obj.d(1))/obj.N_p; 

            % Correction factor:
%             C_M = 0.8; % solid disk gears

            % Gear blank factor:
%             C_R = 1.0; % solid disk gears

            % Basic rack factor:
%             alpha_Pn = obj.alpha_n; % [rad.], Normal pressure angle of basic rack
%             C_B1 = (1.0 + 0.5*(1.2 - obj.h_fP/obj.m_n))*(1.0 - 0.02*(20.0 - alpha_Pn)); % 0.975
%             C_B2 = (1.0 + 0.5*(1.2 - obj.h_fP/obj.m_n))*(1.0 - 0.02*(20.0 - alpha_Pn));

%             C_B = 0.5*(C_B1 + C_B2); % 0.975

            % c' is the Single stiffness:
%             obj.cprime = obj.cprime_th*C_M*C_R*C_B*cosd(obj.beta);

%             % mesh stiffness:
%             obj.c_gamma_alpha = obj.cprime*(0.75*obj.eps_alpha + 0.25);
            
            % Tip relief by running-in:
            C_ay = (1.0/18.0)*(sigma_Hlim/97.0 - 18.45)^2 + 1.5;
            
            K_v = dynamic_factor(obj, n_1, v, rho, u_tmp, F_t*K_gamma*K_A/obj.b, C_ay, f_pbeff, f_falphaeff);
            
            [K_Hbeta,  ~] = obj.face_load_factor();
            
            % Determinant tangential load in a transverse plane:
            F_tH = F_t*K_gamma*K_A*K_v*K_Hbeta;
            
            [K_Halpha, ~] = transv_load_factor(obj, F_tH, f_pb, y_alpha, Z_eps);
            
            % Zone factor: (sec. 6.1)
            Z_H = obj.zone_factor();

            % Single pair tooth contact factor (sec. 6.2)
            [Z_B, Z_D]  = obj.tooth_contact_factor();
            
            % Elasticity factor: (sec. 7)
            Z_E = sqrt(E/(2.0*pi*(1.0 - nu^2)));

            % Helix angle factor: (sec. 9)
            Z_beta = 1.0/sqrt(cosd(obj.beta));

            % Contact stress:
            % Nominal contact stress at pitch point:
            num = F_t*(u_tmp + 1.0);
            den = obj.d(1)*obj.b*u_tmp;
            sigma_H0 = Z_H*Z_E*Z_eps*Z_beta*sqrt(num/den);
            
            % nominal contact stress at pitch point:
            sigma_H1 = Z_B*sigma_H0*sqrt(K_gamma*K_A*K_v*K_Hbeta*K_Halpha); % pinion
            sigma_H2 = Z_D*sigma_H0*sqrt(K_gamma*K_A*K_v*K_Hbeta*K_Halpha); % wheel
            
            [Z_L, Z_v] = lub_vel_factor(sigma_Hlim, nu_40, v);

            % Roughness factor:
            Z_R = obj.rough_factor(R_zh, sigma_Hlim);

            % Work hardening factor:
            Z_W = 1.0;

            % Size factor:
            Z_X = 1.0;

            % Number of load cycles:
            N_L1 = n_1*60.0*L_h; % [-], pinion
            N_L2 = N_L1/u_tmp; % wheel
            
            % Life factor:
            line = 4;
            Z_NT1 = Gear_Set.interp_ZNT(N_L1, line);
            Z_NT2 = Gear_Set.interp_ZNT(N_L2, line);
            
            % Permissible contact stress:
            sigma_HP1 = sigma_Hlim*Z_NT1*Z_L*Z_v*Z_R*Z_W*Z_X/S_Hmin;
            sigma_HP2 = sigma_Hlim*Z_NT2*Z_L*Z_v*Z_R*Z_W*Z_X/S_Hmin;

            % Safety factor for surface durability (against pitting):
            S_H1 = sigma_HP1*S_Hmin/sigma_H1; % pinion/planet
            S_H2 = sigma_HP2*S_Hmin/sigma_H2; % wheel/sun
            
            S_H     = [S_H1 S_H2];
            sigma_H = [sigma_H1 sigma_H2];
        end
        
        function Z_R = rough_factor(obj, R_zh, sigma_Hlim)
            rho_1 = 0.5*obj.d_b(1)*tand(obj.alpha_wt);
            rho_2 = 0.5*obj.d_b(2)*tand(obj.alpha_wt);

            rho_red = (rho_1*rho_2)/(rho_1 + rho_2);
            rho_red = real(rho_red);

        %     R_z = (R_z1 + R_z2)/2.0;
            R_z10 = R_zh*nthroot(10.0/rho_red, 3.0);

            if(sigma_Hlim  < 850.0) % [N/mm^2]
                C_ZR = 0.15;
            elseif((850.0 <= sigma_Hlim) && (sigma_Hlim  < 1200.0))
                C_ZR = 0.32 - sigma_Hlim*2.0e-4;
            else
                C_ZR = 0.08;
            end

            Z_R = power(3.0/R_z10, C_ZR);
        end
        
        function [Z_B, Z_D] = tooth_contact_factor(obj)
            Z_B = NaN;      Z_D = NaN;
            M_1 = tand(obj.alpha_wt)/sqrt((sqrt((obj.d_a(1)/obj.d_b(1))^2 - 1.0) - 2.0*pi/obj.z(1))*(sqrt((obj.d_a(2)/obj.d_b(2))^2 - 1.0) - (obj.eps_alpha - 1.0)*2.0*pi/obj.z(2)));
            M_2 = tand(obj.alpha_wt)/sqrt((sqrt((obj.d_a(2)/obj.d_b(2))^2 - 1.0) - 2.0*pi/obj.z(2))*(sqrt((obj.d_a(1)/obj.d_b(1))^2 - 1.0) - (obj.eps_alpha - 1.0)*2.0*pi/obj.z(1)));

            if((obj.eps_beta == 0.0) && (obj.eps_alpha > 1.0))
                if(M_1 > 1.0)
                    Z_B = M_1;
                else
                    Z_B = 1.0;
                end

                if(M_2 > 1.0)
                    Z_D = M_2;
                else
                    Z_D = 1.0;
                end
            elseif((obj.eps_alpha > 1.0) && (obj.eps_beta >= 1.0))
                Z_B = 1.0;
                Z_D = 1.0;
            elseif((obj.eps_alpha > 1.0) && (obj.eps_beta <  1.0))
                Z_B = M_1 - obj.eps_beta*(M_1 - 1.0);
                Z_D = M_2 - obj.eps_beta*(M_2 - 1.0);
            end
        end
        
        function Z_H = zone_factor(obj)
            num = 2.0*cosd(obj.beta_b)*cosd(obj.alpha_wt);
            den = sind(obj.alpha_wt)*cosd(obj.alpha_t)^2;
            Z_H = sqrt(num/den);
        end
        
        function [K_Halpha, K_Falpha] = transv_load_factor(obj, F_tH, f_pb, y_alpha, Z_eps)
            if(obj.eps_gamma <= 2.0)
                K_Falpha = (0.9 + 0.4*obj.c_gamma_alpha*(f_pb - y_alpha)/(F_tH/mean([obj.b])))*(obj.eps_gamma/2.0);
            else
                K_Falpha = 0.9 + 0.4*sqrt(2.0*(obj.eps_gamma - 1.0)/obj.eps_gamma)*obj.c_gamma_alpha*(f_pb - y_alpha)/(F_tH/mean([obj.b]));
            end
            
            K_Halpha = K_Falpha;
            
            % Transverse load factor (contact stress):
            K_Halpha_lim = obj.eps_gamma/(obj.eps_alpha*Z_eps^2);
            if(K_Halpha > K_Halpha_lim)
                K_Halpha = K_Halpha_lim;
            elseif(K_Halpha < 1.0)
                K_Halpha = 1.0;
            end

            % Transverse load factor (root stress):
            K_Falpha_lim = obj.eps_gamma/(0.25*obj.eps_alpha + 0.75);
            if(K_Falpha > K_Falpha_lim)
                K_Falpha = K_Falpha_lim;
            elseif(K_Falpha < 1.0)
                K_Falpha = 1.0;
            end

        end
        
        function [K_Hbeta, K_Fbeta] = face_load_factor(obj)
            % Face load factor (contact stress): very hard to calculate
            K_Hbeta = 1.15; % according to [2], p. 36
            % h_1 = h_aP + h_fP + k_1*m_n;
            h_1 = abs(obj.d_a(1) - obj.d_f(1))/2.0;        bh1 = obj.b/h_1;
            h_2 = abs(obj.d_a(2) - obj.d_f(2))/2.0;        bh2 = obj.b/h_2;

            bh = min(bh1, bh2);

            if(bh < 3.0)
                bh = 3.0;
            end

            hb = 1.0/bh;

            N_F = 1.0/(1.0 + hb + hb^2);

            % Face load factor (root stress):
            K_Fbeta = power(K_Hbeta, N_F);
        end
        
        function v = pitch_line_vel(obj, n)
            if(strcmp(obj.configuration, "parallel"))
                v = pi*obj.d(1)*n/60.0e3;
            elseif(strcmp(obj.configuration, "planetary"))
                v = (pi*n/60.0e3)*abs(obj.a_w/obj.u - obj.d(1));
%                 v = (pi*n/60.0e3)*abs(obj.a_w - obj.u*obj.d(1));
            else
                error("prog:input", "Configuration [%s] is NOT defined.", obj.configuration);
            end
        end
        
        function K_v = dynamic_factor(obj, n_1, v, rho, u, FtKAb, C_ay, f_pbeff, f_falphaeff)
            
            cond = (v*obj.z(1)/100.0)*sqrt((u^2)/(1.0 + u^2));
            if(cond < 3.0) % [m/s]
                warning("prog:input", "Calculating K_v using method B outside of its useful range. See the end of Sec. 6.3.2 of ISO 6336-1:2006.");
            end
            
            if(strcmp(obj.configuration, "parallel"))
                num = rho*pi*(u*obj.d_m(1)^2)^2;
                den = 8.0*(u^2 + 1.0)*obj.d_b(1)^2;
                m_red = num/den;
%                 (pi/8.0)*(rho*(obj.d_m(1)*u)^2)/(1.0 + u^2)*(obj.d_m(1)/obj.d_b(1))^2;  % solid construction, i.e. 1 - q^4 = 1

                % Resonance running speed:
                n_E1 = 3.0e4*sqrt(obj.c_gamma_alpha/m_red)/(pi*obj.z(1)); % [1/min]

                % Resonance ratio:
                N = n_1/n_E1;
                
                K_v = dyn_factor(FtKAb, N, obj.eps_gamma, obj.cprime, C_ay, f_pbeff, f_falphaeff);
                
            elseif(strcmp(obj.configuration, "planetary"))
                m_sun = (pi/8.0)*rho*(obj.d_m(1)/obj.d_b(1))^4; % sun
                m_pla = (pi/8.0)*rho*(obj.d_m(2)/obj.d_b(2))^4; % planet
                m_red1 = (m_pla*m_sun)/(m_sun + m_pla*obj.N_p);
                m_red2 =  m_pla;

                % Resonance running speed:
                n_E11 = 3.0e4*sqrt(obj.cprime/m_red1)*1.0/(pi*obj.z(1)); % [1/min]
                n_E12 = 3.0e4*sqrt(obj.cprime/m_red2)*1.0/(pi*obj.z(2)); % [1/min]

                % Resonance ratio:
                N_1 =  n_1/n_E11;
                N_2 = (n_1/u)/n_E12;

                K_v1 = dyn_factor(FtKAb, N_1, obj.eps_gamma, obj.cprime, C_ay, f_pbeff, f_falphaeff);
                K_v2 = dyn_factor(FtKAb, N_2, obj.eps_gamma, obj.cprime, C_ay, f_pbeff, f_falphaeff);
                
                K_v = max(K_v1, K_v2);
            else
                error("prog:input", "Configuration [%s] is NOT defined.", obj.configuration);
            end
            
            if(K_v < 1.05) % according to [7], Sec. 7.2.3.2
                K_v = 1.05;
            end
        end
        
        function [S_H, sigma_H] = pitting_factors(obj, P, n_1, varargin)
            if(isempty(varargin))
                mn = obj.m_n;
                bb = obj.b;
            elseif(length(varargin) == 2)
                mn = varargin{1};
                bb = varargin{2};
            elseif(length(varargin{:}) == 1)
                gm = varargin{:};
                mn = obj.m_n*gm;
                bb = obj.b  *gm;
            elseif(length(varargin{:}) == 2) % varargin is a vector with 2 elements
                gm = varargin{:};
                mn = obj.m_n*gm(1);
                bb = obj.b  *gm(2);
            else
                error("prog:input", "Not enough input args. %d ~= 0, 1 or 2.", length(varargin));
            end
            
            aw = obj.a_w*(mn/obj.m_n);
            
            S_Hmin = 1.25;      % [-], Minimum required safety factor for surface durability according to 
            L_h    = 20*365*24; % [h],  Required life
            Q      = 6;         % [-],  ISO accuracy grade
            R_a    = 0.8;       % [um], Maximum arithmetic mean roughness for external gears according to [7], Sec. 7.2.7.2.
            K_A    = 1.25;      % [-], Application factor
            
            mn = Rack.module(mn, "calc", "nearest");

            mod_set = obj.modify_Gear_Set(aw, bb, mn);
            
            [S_H, sigma_H] = mod_set.Pitting_ISO(P, n_1, S_Hmin, L_h, Q, R_a, K_A);
        end
        
        %% Misc.:
        function m_tot = get_mass(obj)
            if(strcmp(obj.configuration, "parallel"))
                m_tot = sum(obj.mass);
            elseif(strcmp(obj.configuration, "planetary"))
                m_tot = sum([1.0, obj.N_p, 1.0].*obj.mass) + obj.carrier.mass;
            end
        end
        
        function J_z = get_J_z(obj)
            % with respect to gear 1.
            if(strcmp(obj.configuration, "parallel"))
                J_z = sum(obj.J_z) + obj.mass(2)*obj.a_w^2;
            elseif(strcmp(obj.configuration, "planetary"))
                J_z = sum(obj.J_z) + obj.N_p*obj.mass(2)*obj.a_w^2 + obj.carrier.J_z;
            end
        end
        
        function aw = find_center_distance(obj, alpha_wt_star)
            aw = abs(obj.z(1) + obj.z(2))*obj.m_n*cosd(obj.alpha_t)/(2.0*cosd(alpha_wt_star)*cosd(obj.beta));
        end
        
        function val = ISO_constraints(obj, varargin)
           if(isempty(varargin))
                aw = obj.a_w;
                bb = obj.b;
                mn = obj.m_n;
            elseif(length(varargin) == 3)
                aw = varargin{1};
                bb = varargin{2};
                mn = varargin{3};
            else
                error("prog:input", "Not enough input args. %d ~= 0 or 3.", length(varargin));
            end
            
            tmp_set = obj.modify_Gear_Set(aw, bb, mn);
            
            val.eps_alp = tmp_set.eps_alpha;
            val.eps_gam = tmp_set.eps_gamma;
            val.d       = tmp_set.d;
        end
        
    end
    
    %% Set methods:
    methods
        function obj = set.a_w(obj, val)
            obj.a_w = val;
        end
        
        function obj = set.cprime(obj, val)
            obj.cprime = val;
        end
        
        function obj = set.shaft(obj, val)
            obj.shaft = val;
        end
        
    end
    
    %% Get methods:
    methods
        function val = get.k_mesh(obj)
            val = obj.c_gamma.*obj.b.*(1.0e6);
        end
        
        function val = get.carrier(obj)
            if(strcmp(obj.configuration, "planetary"))
                val = Carrier(obj.a_w, obj.b);
            else
                error("prog:input", "Only Planetary gear sets have a planet carrier.");
            end
        end
        
        function g = gear(obj, idx)
            g = Gear(obj.m_n, obj.alpha_n, obj.type, obj.z(idx), obj.b, obj.x(idx), obj.beta, obj.k(idx), obj.bore_ratio(idx));
        end
        
        function val = get.u(obj)
            if(strcmp(obj.configuration, "parallel"))
                val = obj.z(2)/obj.z(1);
            elseif(strcmp(obj.configuration, "planetary"))
                val = 1.0 + abs(obj.z(3))/obj.z(1);
            else
                error("prog:input", "Configuration [%s] is NOT defined.", obj.configuration);
            end
        end
        
        function val = get.alpha_wt(obj)
            val = acosd(abs(obj.z(1) + obj.z(2))*(obj.m_n*cosd(obj.alpha_t)/(2.0*obj.a_w*cosd(obj.beta))));
        end
        
        function val = get.eps_alpha(obj)
            B = @(x)(obj.h_fP - x*obj.m_n + obj.rho_fP*(sind(obj.alpha_n) - 1.0));

            d_soi1 = 2.0*sqrt((obj.d(1)/2.0 - B(obj.x(1)))^2 + (B(obj.x(1))/tand(obj.alpha_t))^2);
            d_soi2 = 2.0*sqrt((obj.d(2)/2.0 - B(obj.x(2)))^2 + (B(obj.x(2))/tand(obj.alpha_t))^2);

            % roll angles from the root form diameter to the working pitch point:
            xi_fw1_1 = tand(obj.alpha_wt); % limited by the base diameters
            xi_fw2_1 = tand(obj.alpha_wt);
            
            xi_fw1_2 = tand(obj.alpha_wt) - tand(acosd(obj.d_b(1)/d_soi1)); % limited by the root form diameters
            xi_fw2_2 = tand(obj.alpha_wt) - tand(acosd(obj.d_b(2)/d_soi2));

            xi_fw1_3 = (tand(acosd(obj.d_b(2)/obj.d_a(2))) - tand(obj.alpha_wt))*(obj.z(2)/obj.z(1)); % limited by the the tip diameters of the wheel/pinion
            xi_fw2_3 = (tand(acosd(obj.d_b(1)/obj.d_a(1))) - tand(obj.alpha_wt))*(obj.z(1)/obj.z(2));

            xi_fw1 = min([xi_fw1_1, xi_fw1_2, xi_fw1_3]);
            xi_fw2 = min([xi_fw2_1, xi_fw2_2, xi_fw2_3]);

            % roll angle from the working pitch point to the tip diameter:
            xi_aw1 = xi_fw2*obj.z(2)/obj.z(1);

            % pinion angular pitch:
            tau_1 = 2.0*pi/obj.z(1);

            % Transverse contact ratio:
            val = (xi_fw1 + xi_aw1)/tau_1;

        end
        
        function val = get.eps_beta(obj)
            bb = mean(obj.b);
            val = bb*sind(obj.beta)/(pi*obj.m_n);
        end
        
        function val = get.eps_gamma(obj)
            val = obj.eps_alpha + obj.eps_beta;
        end
        
        function val = get.cprime(obj)
            % Correction factor:
            C_M = 0.8; % solid disk gears

            % Gear blank factor:
            C_R = 1.0; % solid disk gears

            % Basic rack factor:
            alpha_Pn = obj.alpha_n; % [rad.], Normal pressure angle of basic rack
            C_B1 = (1.0 + 0.5*(1.2 - obj.h_fP/obj.m_n))*(1.0 - 0.02*(20.0 - alpha_Pn)); % 0.975
            C_B2 = (1.0 + 0.5*(1.2 - obj.h_fP/obj.m_n))*(1.0 - 0.02*(20.0 - alpha_Pn));

            C_B = 0.5*(C_B1 + C_B2); % 0.975
            
            val = obj.cprime_th*C_M*C_R*C_B*cosd(obj.beta);
        end
        
        function val = get.cprime_th(obj)
            C_1 =  0.04723;      C_2 =  0.15551;      C_3 =  0.25791;
            C_4 = -0.00635;      C_5 = -0.11654;      C_6 = -0.00193;
            C_7 = -0.24188;      C_8 =  0.00529;      C_9 =  0.00182;
            
            % q' is the minimum value for the flexibility of a pair of teeth
            qprime = C_1 + C_2/obj.z_n(1) + C_3/obj.z_n(2) + C_4*obj.x(1) + C_5*(obj.x(1)/obj.z_n(1)) + ...
                C_6*obj.x(2) + C_7*(obj.x(2)/obj.z_n(2)) + C_8*obj.x(1)^2 + C_9*obj.x(1)^2; % [mm-um/N]
            
            % c'_th is the theoretical single stiffness:
            val = 1.0/qprime;

        end
        
        function val = get.c_gamma(obj)
            val = obj.c_gamma_alpha + obj.c_gamma_beta;
        end
        
        function val = get.c_gamma_alpha(obj)
            val = obj.cprime*(0.75*obj.eps_alpha + 0.25);
        end
        
        function val = get.c_gamma_beta(obj)
            val = 0.85*obj.c_gamma_alpha;
        end
        
    end
    
end

function Z_eps = contact_ratio_factor(beta, eps_alpha, eps_beta)

    if(beta == 0.0)
        Z_eps = sqrt((4.0 - eps_alpha)/3.0);
    else
        if(eps_beta < 1.0)
            Z_eps = sqrt((1.0 - eps_beta)*(4.0 - eps_alpha)/3.0 + eps_beta/eps_alpha);
        else
            Z_eps = sqrt(1.0/eps_alpha);
        end
    end
end

function K_v = dyn_factor(FtKAb, N, eps_gamma, cp, C_a, f_pbeff, f_falphaeff)

    if(FtKAb < 100) % [N/mm]
        N_S = 0.5 + 0.35*sqrt(FtKAb/100.0);
    else
        N_S = 0.85;
    end

    C_v1 = 0.32;        C_v5 = 0.47;

    if((1.0 < eps_gamma) && (eps_gamma <= 2.0))
        C_v2 = 0.34;        C_v3 = 0.23;
        C_v4 = 0.90;        C_v6 = 0.47;
    elseif(eps_gamma > 2.0)
        C_v2 =  0.57 /(eps_gamma - 0.3);
        C_v3 =  0.096/(eps_gamma - 1.56);
        C_v4 = (0.57 - 0.05*eps_gamma)/(eps_gamma - 1.44);
        C_v6 =  0.12/(eps_gamma - 1.74);
    else
        C_v2 = NaN;     C_v3 = NaN;
        C_v4 = NaN;     C_v6 = NaN;
    end
    
    if((1.0 < eps_gamma) && (eps_gamma <= 1.5))
        C_v7 = 0.75;
    elseif((1.5 < eps_gamma) && (eps_gamma <= 2.5))
        C_v7 = 0.125*sind(pi*(eps_gamma - 2.0)) + 0.875;
    elseif(eps_gamma > 2.5)
        C_v7 = 1.0;
    end

    B_p = cp*f_pbeff/FtKAb;
    B_f = cp*f_falphaeff/FtKAb;
    B_k = abs(1.0 - cp*C_a/FtKAb);

    % Dynamic factor:
    K_v = NaN;
    if(N <= N_S)
        K = C_v1*B_p + C_v2*B_f + C_v3*B_k; 
        K_v = N*K + 1.0;
    elseif((N_S < N) && (N <= 1.15))
        K_v = C_v1*B_p + C_v2*B_f + C_v4*B_k + 1.0;
    elseif((1.15 < N) && (N < 1.5))
        K_vN115 = C_v1*B_p + C_v2*B_f + C_v4*B_k + 1.0;
        K_vN15  = C_v5*B_p + C_v6*B_f + C_v7;
        K_v = K_vN15 + (K_vN115 - K_vN15)*(1.5 - N)/0.35;
    elseif(N >= 1.5)
        K_v = C_v5*B_p + C_v6*B_f + C_v7;
    end
end

function [Z_L, Z_v] = lub_vel_factor(sigma_Hlim, v_40, v)
    
    if(sigma_Hlim  < 850.0) % [N/mm^2]
        C_ZL = 0.83;
    elseif((850.0 <= sigma_Hlim) && (sigma_Hlim  < 1200.0))
        C_ZL = sigma_Hlim/4375.0 + 0.6357;
    else
        C_ZL = 0.91;
    end

    % Lubricant factor:
    Z_L = C_ZL + 4.0*(1.0 - C_ZL)/(1.2 + 134.0/v_40)^2;

    % Velocity factor:
    C_Zv = C_ZL + 0.02;
    Z_v = C_Zv + 2.0*(1.0 - C_Zv)/sqrt(0.8 + 32.0/v);

end
