classdef Gear < Rack
    %GEAR This class implements the geometric concepts and parameters for
    % cylindrical gears with involute helicoid tooth flanks. It also
    % implements the concepts and parameters for cylindrical gear pairs
    % with parallel axes and a constant gear ratio.
    % References: 
    % [1] ISO 21771:2007 Gears -- Cylindrical involute gears and gear pairs
    % -- Concepts and geometry
    % [2] ISO 1328-1:1995 Cylindrical gears -- ISO system of accuracy --
    % Part 1: Definitions and allowable values of deviations relevant to
    % corresponding flanks of gear teeth
    % [3] ISO 1122-1: Vocabulary of gear terms -- Part 1: Definitions
    % related to geometry.
    % 
    % Definitions:
    % - Driving gear: that gear of a gear pair which turns the other.
    % - Driven gear: that gear of a gear pair which is turned by the other.
    % - Pinion: that gear of a pair which has the smaller number of teeth.
    % - Wheel: that gear of a pair which has the larger number of teeth.
    % - Gear ratio: quotient of the number of teeth of the wheel divided by
    % the number of teeth of the pinion.
    % - Transmission ratio: quotient of the angular speed of the first
    % driving gear divided by the angular speed of the last driven gear of
    % a gear train.
    % - (Right/Left)-handed teeth: teeth whose sucessive transverse 
    % profiles show (clockwise/anti-clockwise) displacement with increasing
    % distance from an observer looking along the straight line generators
    % of the reference surface.
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
    
    properties
        z         (1, :) {mustBeInteger, mustBeFinite}                    = 13;  % [-],    Number of teeth
        x         (1, :) {mustBeNumeric, mustBeFinite}                    = 0.0; % [-],    Profile shift coefficient
        beta      (1, :) {mustBeNumeric, mustBeFinite, mustBeNonnegative} = 0.0; % [deg.], Helix angle (at reference cylinder)
        k         (1, :) {mustBeNumeric, mustBeFinite}                    = 0.0; % [-],    Tip alteration coefficient
        bore_ratio(1, :) {mustBeNumeric, mustBeFinite, mustBePositive}    = 0.5; % [-],    Ratio btw. bore and reference diameters
        Q         (1, 1) {mustBeNumeric, mustBeFinite, mustBePositive}    = 6.0; % [-],    ISO accuracy grade
        R_a       (1, 1) {mustBeNumeric, mustBeFinite, mustBePositive}    = 1.0; % [um],   Arithmetic mean roughness
        material  (1, :) Material;
    end
    
    properties(Access = public)
        b         (1, :) {mustBeNumeric, mustBeFinite, mustBePositive}    = 13;  % [mm],   Face width
    end
    
    properties(Access = private)
        m_int;
        d_int;
        b_int;
    end
    
    properties(Dependent)
        m_n;     % [mm],     Normal module
        alpha_n; % [deg.],   Pressure angle (at reference cylinder)
        m_t;     % [mm],     Transverse module
        alpha_t; % [rad.],   Transverse pressure angle
        beta_b;  % [deg.],   Base helix angle
        p_t;     % [mm],     Transverse pitch
        p_bt;    % [mm],     Transverse base pitch
        p_et;    % [mm],     Transverse base pitch on the path of contact
        h;       % [mm],     Tooth depth
        d;       % [mm],     Reference diameter
        d_a;     % [mm],     Tip diameter
        d_b;     % [mm],     Base diameter
        d_f;     % [mm],     Root diameter
        d_m;     % [mm],     Mean tooth diameter
        d_w;     % [mm],     Working pitch diameter
        d_bore;  % [mm],     Bore diameter
        z_n;     % [-],      Virtual number of teeth
        V;       % [m^3],    Volume
        mass;    % [kg],     Mass
        J_x;     % [kg-m^2], Mass moment of inertia (rot. axis)
        J_y;     % [kg-m^2], Mass moment of inertia
        J_z;     % [kg-m^2], Mass moment of inertia
        f_pt;    % [um],     Single pitch deviation according to Sec. 6.1 of ISO 1328-1 [2]
        F_p;     % [um],     Total cumulative pitch deviation according to Sec. 6.3 of ISO 1328-1 [2]
        F_alpha; % [um],     Total profile deviation according to Section 6.4 of ISO 1328-1 [2]
        f_beta;  % [um],     Total helix deviation according to Section 6.5 of ISO 1328-1 [2]
        f_falpha;% [um],     Profile form deviation according to App. B.2.1 of ISO 1328-1 [2]
        f_Halpha;% [um],     Profile slope deviation according to App. B.2.2 of ISO 1328-1 [2]
        f_fbeta; % [um],     Helix form deviation according to App. B.2.3 of ISO 1328-1 [2]
        f_Hbeta; % [um],     Helix slope deviation according to App. B.2.3 of ISO 1328-1 [2]
        d_Ff;    % [mm],     Root form diameter
        R_z;     % [um],     Mean peak-to-valley surface roughness
    end
    
    methods
        function obj = Gear(varargin)
            default = {'m_n'        ,  8.0,  ...
                       'alpha_n'    ,  20.0,  ...
                       'type'       , 'D',  ...
                       'z'          ,  17,   ...
                       'b'          , 100.0, ...
                       'x'          ,   0.145,  ...
                       'beta'       ,  15.8,  ...
                       'k'          ,   0.0,  ...
                       'bore_ratio' ,   0.5,  ...
                       'Q'          ,   5.0,  ...
                       'R_a'        ,   1.0, ...
                       'material'   ,   Material()};
            
            default = process_varargin(default, varargin);
            
            obj@Rack('type'   , default.type, ...
                     'm'      , default.m_n, ...
                     'alpha_P', default.alpha_n);

            obj.z          = default.z;
            obj.b          = default.b;
            obj.x          = default.x;
            obj.beta       = default.beta;
            obj.k          = default.k;
            obj.bore_ratio = default.bore_ratio;
            obj.Q          = default.Q;
            obj.R_a        = default.R_a;
            obj.material   = default. material;
            
            range_b = [  0.0,   4.0, 10.0, 20.0, 40.0, 80.0, 160.0, 250.0, ...
                       400.0, 650.0,  1.0e3];
            range_d = [0.0,   5.0,  20.0,  50.0, 125.0,  280.0, 560.0, 1.0e3, ...
                       1.6e3, 2.5e3, 4.0e3, 6.0e3, 8.0e3, 10.0e3];
            range_m = [ 0.0,  0.5, 2.0, 3.5, 6.0, 10.0, 16.0, 25.0,  ...
                       40.0, 70.0];

            idx_b = find(range_b > obj.b, 1, 'first');
            obj.b_int = geomean(range_b(idx_b-1:idx_b));
            idx_m = find(range_m > obj.m_n, 1, 'first');
            obj.m_int = geomean(range_m(idx_m-1:idx_m));
            
            for idx = 1:length(obj.d)
                dd  = obj.d(idx);
                jdx = find(range_d > dd, 1, 'first');
                obj.d_int(idx) = geomean(range_d(jdx-1:jdx));
            end
        end
        
        function tab = disp(obj)
            tab_str = {"Normal module",                                "m_n",     "mm",     obj.m_n;
                       "Pressure angle",                               "alpha_n", "deg.",   obj.alpha_n;
                       "Number of teeth",                              "z",       "-",      obj.z;
                       "Face width",                                   "b",       "mm",     obj.b;
                       "Profile shift coefficient",                    "x",       "-",      obj.x;
                       "Helix angle",                                  "beta",    "deg.",   obj.beta;
                       "Tip alteration coefficient",                   "k",       "-",      obj.k;
                       "Transverse module",                            "m_t",     "mm",     obj.m_t;
                       "Transverse pressure angle",                    "alpha_t", "deg.",   obj.alpha_t;
                       "Base helix angle",                             "beta_b",  "deg.",   obj.beta_b;
                       "Normal pitch",                                 "p_n",     "mm",     obj.p;
                       "Transverse pitch",                             "p_t",     "mm",     obj.p_t;
                       "Transverse base pitch",                        "p_bt",    "mm",     obj.p_bt;
                       "Transverse base pitch on the path of contact", "p_et",    "mm",     obj.p_et;
                       "Reference diameter",                           "d",       "mm",     obj.d;
                       "Tip diameter",                                 "d_a",     "mm",     obj.d_a;
                       "Base diameter",                                "d_b",     "mm",     obj.d_b;
                       "Root diameter",                                "d_f",     "mm",     obj.d_f;
                       "Mean tooth diameter",                          "d_m",     "mm",     obj.d_m;
                       "Mass",                                         "m",       "kg",     obj.mass;
                       "Mass moment of inertia (x axis, rot.)",        "J_x",     "kg-m^2", obj.J_x;
                       "Mass moment of inertia (y axis)",              "J_y",     "kg-m^2", obj.J_y;
                       "Mass moment of inertia (z axis)",              "J_z",     "kg-m^2", obj.J_z;
                       };

            Parameter =          tab_str(:,1);
            Symbol    =          tab_str(:,2);
            Unit      =          tab_str(:,3);
            Value     = cell2mat(tab_str(:,4));
%             Value     = tab_str(:,4);

            tab = table(Parameter, Symbol, Value, Unit, ...
                        'variableNames', ["Parameter", "Symbol", "Value", "Unit"]);
            
            if(nargout == 0)
                disp(tab);
                clear tab;
            end
        end
        
        function h = plot(obj, varargin)
            if(nargin == 1)
                C = zeros(2, 1);
                plot_prop = {'lineStyle', '-' , 'lineWidth', 2.0, 'color', [1.0 0.0 0.0]};
            elseif(nargin > 1)
                C = varargin{1};
                plot_prop = varargin(2:end);
            else
                error('prog:input', 'Too many variables.');
            end
            
            [X, Y] = obj.reference_circle(C);
            
            axis equal;
            box on;
            h = plot(X, Y, plot_prop{:});
            
        end
        
        function [X, Y, Z] = plot3(obj, varargin)
            % adapted from: (access on 20/11/2019)
            % https://www.mathworks.com/matlabcentral/answers/62894-trying-to-plot-a-3d-closed-cylinder
            
            if(nargin == 1)
                C = zeros(2, 1);
                plot_prop = {'edgeColor', 'none', 'lineStyle', '-' , 'faceColor', [1.0 0.0 0.0]};
            elseif(nargin > 1)
                C = varargin{1};
                plot_prop = varargin(2:end);
            else
                error('prog:input', 'Too many variables.');
            end
            
            if(obj.z > 0)
                [X_tmp, Y_tmp, Z_tmp] = cylinder(obj.d/2, obj.z);
                X =  X_tmp + C(1);
                Y =  Z_tmp.*obj.b;
                Z = -Y_tmp - C(2);

                surf(X, Y, Z, plot_prop{:});
                hold on;
                fill3(X(1,:), Y(1,:), Z(1,:), plot_prop{end});
                fill3(X(2,:), Y(2,:), Z(2,:), plot_prop{end});
            else
                % External cylinder
                [X_tmp, Y_tmp, Z_tmp] = cylinder(obj.d_bore/2, abs(obj.z));
                X =  X_tmp + C(1);
                Y =  Z_tmp.*obj.b;
                Z = -Y_tmp - C(2);

                % Internal cylinder
                [X2_tmp, Y2_tmp, Z2_tmp] = cylinder(obj.d/2, abs(obj.z));
                X2 =  X2_tmp + C(1);
                Y2 =  Z2_tmp.*obj.b;
                Z2 = -Y2_tmp - C(2);

                [x1_tmp, y1_tmp, z1_tmp] = fill_ring(obj.d_bore/2, obj.d/2, abs(obj.z));
                x1 =  x1_tmp + C(1);
                y1 =  z1_tmp;
                z1 = -y1_tmp - C(2);
                
                surf(X, Y, Z, plot_prop{:});
                hold on;
                
                surf(X2, Y2,         Z2, plot_prop{:});
                surf(x1, y1        , z1, plot_prop{:});
                surf(x1, y1 + obj.b, z1, plot_prop{:});
            end
            
            hold off;
            box on;
            axis equal;
            axis ij;
            
            if(nargout == 0)
                clear X Y Z;
            end
        end
        
        function h = rectangle(obj, varargin)
            if(nargin == 1)
                C = zeros(2,1);
                plot_prop = {[1.0 0.0 0.0], 'edgeColor', 'k', 'lineStyle', '-' , 'faceColor', [1.0 0.0 0.0]};
            elseif(nargin > 1)
                C = varargin{1};
                plot_prop = varargin(2:end);
            else
                error('prog:input', 'Too many variables.');
            end
            
            X = 0.5.*obj.b.*[1 -1 -1  1] + C(1);
            Y = 0.5.*obj.d.*[1  1 -1 -1] + C(2);
            h = fill(X, Y, plot_prop{:});
            
            axis equal;
            box on;
            
        end
    end
    
    %% Calculations:
    methods
        function f_pt = single_pitch_tol(obj, Q)
        %SINGLE_PITCH_TOL Single pitch tolerance, according to ISO
        % 1328-1:1995.

            d_list = [5.0 20.0 50.0 125.0 280.0 560.0 1000.0 1600.0 2500.0 4000.0 6000.0 8000.0 10000.0];
            m_list = [0.5 2.0 3.5 6.0 10.0 16.0 25.0 40.0 70.0];
            
            idx_m = find(m_list > obj.m_n, 1, 'first');
            m_int = geomean(m_list(idx_m-1:idx_m));
            
            f_pt = zeros(size(obj.d));
            
            for idx = 1:length(obj.d)
                idx_d = find(d_list > obj.d(idx)  , 1, 'first');

                d_int = geomean(d_list(idx_d-1:idx_d));

                f_pt(idx) = (0.3.*(m_int + 0.4.*sqrt(d_int)) + 4.0).*power(2.0, (Q - 5.0)/2.0);
                f_pt(idx) = round_ISO(f_pt(idx));
            end

        end
        
        function f_falpha = profile_form_tol(obj, Q)
        %PROFILE_FORM Profile form tolerance, according to ISO 1328-1:1995.

            d_list = [5.0 20.0 50.0 125.0 280.0 560.0 1000.0 1600.0 2500.0 4000.0 6000.0 8000.0 10000.0];
            m_list = [0.5 2.0 3.5 6.0 10.0 16.0 25.0 40.0 70.0];

            idx_m = find(m_list > obj.m_n, 1, 'first');
            m_int = geomean(m_list(idx_m-1:idx_m));
            
            f_falpha = zeros(size(obj.d));

            for idx = 1:length(obj.d)
                idx_d = find(d_list > obj.d(idx), 1, 'first');

                d_int = geomean(d_list(idx_d-1:idx_d));

                f_falpha(idx) = (2.5.*sqrt(m_int) + 0.17.*sqrt(d_int) + 0.5).*power(2.0, (Q - 5.0)/2.0);
                f_falpha(idx) = round_ISO(f_falpha(idx));
            end
            
        end
        
        function obj = work_pitch_diam(obj, alpha_wt)
            %WORK_PITCH_DIAM Working pitch diameter, [mm]
            obj.d_w = obj.d_b./cosd(alpha_wt);
        end
        
        function obj = get_mass(obj, rho)
            %GET_MASS Update the gear's mass with a user defined density.
            if(nargin == 0)
                rho = obj.material.rho*1.0e9;
            end
            obj.mass = rho.*obj.V;
        end
        
        function [x, y] = reference_circle(obj, varargin)
            %REFERENCE_CIRCLE Returns the points to draw the reference
            % circle of a gear.
            t = linspace(0.0, 2.0*pi, 1001);
            
            R = obj.d/2.0;
            
            if(nargin == 1)
                C = zeros(2, 1);
            elseif(nargin == 2)
                C = varargin{1};
            else
                error('prog:input', 'Too many variables.');
            end
    
            x = R*cos(t) + C(1);
            y = R*sin(t) + C(2);

            x = [x, C(1)];
            y = [y, C(2)];
            
        end
        
        function val = round_ISO(obj, x)
            % Credits for the original version of this method go to:
            % E. M. F. Donéstevez, “Python library for design of spur and 
            % helical gears transmissions.” Zenodo, 09-Feb-2020, 
            % doi: 10.5281/ZENODO.3660527.
            %
            % It was modified to account for the case where x is an array 
            % and to remove the second argument.
            %

            x = x.*power(2.0, (obj.Q - 5.0)/2.0);
            val = zeros(size(x));
            
            for idx = 1:length(x)
                xx = x(idx);
                if(xx >= 10.0)
                    val(idx) = round(xx);
                elseif((5.0 <= xx) && (xx <= 10.0))
                    y = mod(mod(xx ,1), 1);
                    if((mod(xx, 1) <= 0.25) || ((0.5 <= y) && (y <= 0.75)))
                        val(idx) = floor(2.0.*xx)/2.0;
                    else
                        val(idx) = ceil(2.0.*xx)/2.0;
                    end
                else
                    val(idx) = round(xx, 1);
                end
            end
            
        end
    end
    
    %% Set methods:
    methods
        function obj = set.m_n(obj, val)
            obj.m = Rack.module(val, 'calc', 'nearest');
        end
        
        function obj = set.b(obj, val)
            obj.b = val;
        end
    end
    
    %% Get methods:
    methods
        function val = get.m_n(obj)
            %M_N Normal module, [mm]
            val = obj.m;
        end
        
        function val = get.alpha_n(obj)
            %ALPHA_N Pressure angle (at reference cylinder), [deg.]
            val = obj.alpha_P;
        end
        
        function val = get.m_t(obj)
            %M_T Transverse module, [mm]
            val = obj.m_n./cosd(obj.beta);
        end
        
        function val = get.alpha_t(obj)
            %ALPHA_T Transverse pressure angle, [rad.]
            val = atand(tand(obj.alpha_n)./cosd(obj.beta));
        end
        
        function val = get.beta_b(obj)
            % [deg.],   Base helix angle
            val = atand(tand(obj.beta).*cosd(obj.alpha_t));
        end
        
        function val = get.p_t(obj)
            % [mm],     Transverse pitch
            val = pi*obj.m_t;
        end
        
        function val = get.p_bt(obj)
            % [mm],     Transverse base pitch
            val = obj.p_t.*cosd(obj.alpha_t);
        end
        
        function val = get.p_et(obj)
            % [mm],     Transverse base pitch on the path of contact
            val = obj.p_bt;
        end
        
        function val = get.h(obj)
            % [mm],     Tooth depth
            val = obj.h_aP + obj.k.*obj.m_n + obj.h_fP;
        end
        
        function val = get.d(obj)
            % [mm],     Reference diameter
            val = abs(obj.z).*obj.m_t;
        end
        
        function val = get.d_a(obj)
            % [mm],     Tip diameter
            val = obj.d + 2.0*sign(obj.z).*(obj.x.*obj.m_n + obj.h_aP + obj.k.*obj.m_n);
        end
        
        function val = get.d_b(obj)
            % [mm],     Base diameter
            val = obj.d.*cosd(obj.alpha_t);
        end
        
        function val = get.d_f(obj)
            % [mm],     Root diameter
            val = obj.d - 2.0*sign(obj.z).*(obj.h_fP - obj.x.*obj.m_n);
        end
        
        function val = get.d_m(obj)
            % [mm],     Mean tooth diameter
            val = (obj.d_a + obj.d_f)/2.0;
        end
        
        function val = get.d_bore(obj)
            % [mm],     Bore diameter
            val = obj.bore_ratio.*obj.d;
        end
        
        function val = get.z_n(obj)
            % [-],      Virtual number of teeth
            val = obj.z/(cosd(obj.beta).*cosd(obj.beta_b).^2);
        end
        
        function val = get.V(obj)
            % [mm^3],   Volume
            val = zeros(size(obj.z));
            
            for idx = 1:length(obj.z)
                if(obj.z(idx) > 0)
                    r_out = (obj.d_a(idx) + obj.d_f(idx))/4.0;
                    r_in  = obj.d_bore(idx)/2.0;
                else
                    r_in = (obj.d_a(idx) + obj.d_f(idx))/4.0;
                    r_out = obj.d_bore(idx)/2.0;
                end

                val(idx) = obj.b.*pi.*(r_out.^2 - r_in.^2);
            end
            val = val*1.0e-9;
        end
        
        function val = get.mass(obj)
            % [kg],     Mass
            rho = obj.material.rho*1.0e9;
            val = rho.*obj.V;
        end
        
        function val = get.J_x(obj)
            % [kg-m^2], Mass moment of inertia (rot. axis)
            if(obj.z > 0)
                r_out = (obj.d_a + obj.d_f)/4.0;
                r_in  = obj.d_bore/2.0;
            else
                r_in = (obj.d_a + obj.d_f)/4.0;
                r_out = obj.d_bore/2.0;
            end
            val = (obj.mass/2.0).*(r_out.^2 + r_in.^2)*1.0e-6;
        end
        
        function val = get.J_y(obj)
            % [kg-m^2], Mass moment of inertia
            if(obj.z > 0)
                r_out = (obj.d_a + obj.d_f)/4.0;
                r_in  = obj.d_bore/2.0;
            else
                r_in = (obj.d_a + obj.d_f)/4.0;
                r_out = obj.d_bore/2.0;
            end
            val = (obj.mass/12.0).*(3.0.*(r_out.^2 + r_in.^2) + obj.b.^2)*1.0e-6;
        end
        
        function val = get.J_z(obj)
            % [kg-m^2], Mass moment of inertia
            val = obj.J_y;
        end
        
        function val = get.f_pt(obj)
            val = 0.3*(obj.m_int + 0.4*sqrt(obj.d_int)) + 4.0;
            val = obj.round_ISO(val);
        end
        
        function val = get.F_p(obj)
            val = 0.3*obj.m_int + 1.25*sqrt(obj.d_int) + 7.0;
            val = obj.round_ISO(val);
        end
        
        function val = get.F_alpha(obj)
            val = 3.2*sqrt(obj.m_int) + 0.22*sqrt(obj.d_int) + 0.7;
            val = obj.round_ISO(val);
        end
        
        function val = get.f_beta(obj)
            val = 0.1*sqrt(obj.d_int) + 0.63*sqrt(obj.b_int) + 4.2;
            val = obj.round_ISO(val);
        end
        
        function val = get.f_falpha(obj)
            val = 2.5*sqrt(obj.m_int) + 0.17*sqrt(obj.d_int) + 0.5;
            val = obj.round_ISO(val);
        end
        
        function val = get.f_Halpha(obj)
            val = 2.0*sqrt(obj.m_int) + 0.14*sqrt(obj.d_int) + 0.5;
            val = obj.round_ISO(val);
        end
        
        function val = get.f_fbeta(obj)
            val = (0.07*sqrt(obj.d_int) + 0.45*sqrt(obj.b_int) + 3.0);
            val = obj.round_ISO(val);
        end
        
        function val = get.f_Hbeta(obj)
            val = obj.f_fbeta;
        end
        
        function val = get.d_Ff(obj)
            % ISO 21771, Sec. 7.6, Eq. (128)
            B = (obj.h_fP - obj.x.*obj.m_n + obj.rho_fP.*(sind(obj.alpha_n) - 1.0));
            val = sqrt((obj.d.*sind(obj.alpha_t) - 2.0.*B./sind(obj.alpha_t)).^2 + obj.d_b.^2);
        end
        
        function val = get.R_z(obj)
            val = 6.0.*obj.R_a;
        end
    end
    
end

function x_r = round_ISO(x)
% rounding according to ISO 1328-1, Sec. 5.2.3:
    
    if(x < 5.0) % [um]
        x_r = round_near(x, 0.1);
    elseif((5.0 <= x) && (x <= 10.0))
        x_r = round_near(x, 0.5);
    else
        x_r = round(x);
    end
end

function x_r = round_near(x, y)
%ROUND_NEAR rounds the number x to the nearest y. Based on:
% https://www.mathworks.com/matlabcentral/answers/14495-rounding-elements-in-array-to-nearest-0-25
% (Accessed on 08/10/2019)
%

    % decimal part:
    dp = x - floor(x);
    
    up =  ceil(dp/y)*y;
    dw = floor(dp/y)*y;
    
    d_up = up - dp;
    d_dw = dp - dw;
    
    dx = zeros(size(x));
    
    for idx = 1:length(x)
        if(d_up(idx) < d_dw(idx))
            dx(idx) = up(idx);
        else
            dx(idx) = dw(idx);
        end
    end
    
    x_r = floor(x) + dx;

end

function [x, y, z] = fill_ring(r_out, r_in, n)
    % https://www.mathworks.com/matlabcentral/answers/178464-how-to-make-circular-ring
    
    t = linspace(0.0, 2.0*pi, n + 1);
    r_outer = ones(size(t))*r_out;
    r_inner = ones(size(t))*r_in;
    r = [r_outer; r_inner];
    t = [t; t];
    x = r.*cos(t);
    y = r.*sin(t);
    z = zeros(size(x));
    
    if(nargout == 0)
        surf(x, y, z);
        clear x y z;
    end
    
end
