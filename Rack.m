classdef Rack
    %RACK This class implements the characteristics of the standard basic 
    % rack tooth profile for cylindrical involute gears (external or 
    % internal) for general and heavy engineering.
    %
    % References:
    % [1] ISO 53:1998 Cylindrical gears for general and heavy engineering
    % -- Standard basic rack tooth profile
    % [2] ISO 54:1996 Cylindrical gears for general engineering and for
    % heavy engineering -- Modules
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
    
    properties(Access = protected)
        m       (1, 1) double {mustBePositive, mustBeInRange(m, [0.01 50.0])} = 1.0;  % [mm],   Module
    end
    
    properties(SetAccess = private)
        type    (1, :) string {mustBeMember(type, ["A", "B", "C", "D"])}     = "A";  % [-],    Type of basic rack tooth profile
        alpha_P (1, 1) double {mustBeNonnegative}                            = 20.0; % [deg.], Pressure angle
        U_FP    (1, 1) double {mustBeNumeric}                                = 0.0;  % [mm],   Size of undercut
        alpha_FP(1, 1) double {mustBeNumeric}                                = 0.0;  % [deg.], Angle of undercut
    end
    
    properties(Access = protected, Dependent)
        c_P;      % [mm], Bottom clearance
        e_P;      % [mm], Spacewidth
        h_aP;     % [mm], Addendum
        h_fP;     % [mm], Dedendum
        h_FfP;    % [mm], Straight portion of the dedendum
        h_P;      % [mm], Tooth depth
        h_wP;     % [mm], Common depth of rack and tooth
        p;        % [mm], Pitch
        s_P;      % [mm], Tooth thickness
        rho_fP;   % [mm], Fillet radius
    end
    
    methods
        function obj = Rack(varargin)
            %RACK creates a Rack object with type, module and pressure
            % angle alpha_P.
            default = {'type'   , 'A', ...
                       'm'      , 1.0, ...
                       'alpha_P', 20.0};
            
            default = scaling_factor.process_varargin(default, varargin);
            
            % fit the module to [2]:
            m = Rack.module(default.m, "round_0125", "nearest");
            
            obj.type     = default.type;
            obj.m        = m;
            obj.alpha_P  = default.alpha_P;
            obj.U_FP     = 0.0;
            obj.alpha_FP = 0.0;
        end
        
        function tab = disp(obj)
            tab_str = {"Module",                           "m",         "mm",   obj.m;
                       "Pressure angle",                   "alpha_P",   "deg.", obj.alpha_P;
                       "Bottom clearance",                 "c_P",       "mm",   obj.c_P;
                       "Space width",                      "e_P",       "mm",   obj.e_P;
                       "Addendum",                         "h_aP",      "mm",   obj.h_aP;
                       "Dedendum",                         "h_fP",      "mm",   obj.h_fP;
                       "Straight portion of the dedendum", "h_FfP",     "mm",   obj.h_FfP;
                       "Tooth depth",                      "h_P",       "mm",   obj.h_P;
                       "Common depth of rack and tooth",   "h_wP",      "mm",   obj.h_wP;
                       "Pitch",                            "p",         "mm",   obj.p;
                       "Tooth thickness",                  "s_P",       "mm",   obj.s_P;
                       "Size of undercut",                 "U_FP",      "mm",   obj.U_FP;
                       "Angle of undercut",                "alpha_FP",  "deg.", obj.alpha_FP;
                       "Fillet radius",                    "rho_fP",    "mm",   obj.rho_fP;
                       "Maximum fillet radius",            "rho_fPmax", "mm",   obj.max_fillet_radius;
                        };
            
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
        
    end
    
    %% Calculations:
    methods(Static)
        function m = module(m_x, option, round_opt)
        %MODULE Returns the values of normal modules for straight and helical
        % gears according to ISO 54:1996. According to this standard, preference
        % should be given to the use of the normal modules as given in series I and
        % the module 6.5 in series II should be avoided.
        % The module is defined as the quotient between:
        % - the pitch, expressed in millimetres, to the number pi, or;
        % - the reference diameter, expressed in millimetres, by the number of teeth.
        %

            m_1 = [1.0 1.25 1.5 2.0 2.5 3.0 4.0 5.0 6.0 8.0 10.0 12.0 16.0 20.0 25.0 32.0 40.0 50.0];
            m_2 = [1.125 1.375 1.75 2.25 2.75 3.5 4.5 5.5 6.5 7.0 9.0 11.0 14.0 18.0 22.0 28.0 36.0 45.0];
%             mini = [0.125,0.150,0.200,0.250,0.300,0.375,0.400,0.500,0.750,0.800];
            switch option
                case "show"
                    tab = table(m_1', m_2');
                    tab.Properties.VariableNames = ["Serie_I", "Serie_II"];
                    tab.Properties.VariableUnits = ["mm", "mm"];
                    disp(tab);
                    fprintf("Obs.: According to ISO 54:1996, preference should be \ngiven to the use of the normal modules as given in \nseries I and the module 6.5 in series II should be avoided.\n");
                case "calc"
                    m_2(9) = 21.0; % removing the value 6.5 which should be avoided according to ISO 54:1996
                    m_2(end + 1) = 30.0;
                    x = sort([m_1 m_2]);
                case "calc_1"
                    x = m_1;
                case "calc_2"
                    m_2(9) = 21.0; % removing the value 6.5 which should be avoided according to ISO 54:1996
                    x = m_2;
                case 'round_0125'
                    
                otherwise
                    error("prog:input", "Option [%s] is NOT valid.", upper(option));
            end
            
            m = zeros(size(m_x));
            
            if(strcmp(option, 'round_0125'))
                m = round(m_x*8)/8;
            elseif(contains(option, 'calc', 'IgnoreCase', true))
                idx = interp1(x, 1:length(x), m_x, round_opt);
                
                for jdx = 1:length(m_x)
                    if(isnan(idx(jdx)))
                        if(m_x(jdx) < x(1))
                            m(jdx) = x(1);
                        elseif(m_x(jdx) > x(end))
                            m(jdx) = x(end);
                        end
                    else
                        m(jdx) = x(idx(jdx));
                    end
                end
                
            end
            
            if(nargout == 0)
                clear m;
            end
        end
        
    end
    
    methods
        function rho_fP_max = max_fillet_radius(obj)
            %MAX_FILLET_RADIUS returns the maximum fillet radius of the
            %basic rack according to ISO 53:1998, Sec. 5.9.
            if(obj.alpha_P ~= 20.0)
                error("prog:input", "The pressure angle is not 20.0 [deg.]");
            else
                if((obj.c_P <= 0.295*obj.m) && (obj.h_FfP == obj.m))
                    rho_fP_max = obj.c_P./(1.0 - sind(obj.alpha_P));
                elseif((0.295*obj.m < obj.c_P) && (obj.c_P <= 0.396*obj.m))
                    rho_fP_max = (pi.*obj.m./4.0 - obj.h_fP.*tand(obj.alpha_P))./tand((90.0 - obj.alpha_P)./2.0);
                end
            end
        end
    end
    
    %% Get methods:
    methods
        function c_P = get.c_P(obj)
            %C_P Bottom clearance, [mm]
            switch obj.type
                case "A"
                    c_P = 0.25*obj.m;
                case "B"
                    c_P = 0.25*obj.m;
                case "C"
                    c_P = 0.25*obj.m;
                case "D"
                    c_P = 0.40*obj.m;
                otherwise
                    error("prog:input", "Type [%s] is NOT valid.", obj.type);
            end
        end
        
        function e_P = get.e_P(obj)
            %E_P Spacewidth, [mm]
            e_P = pi*obj.m/2.0;
        end
        
        function h_aP = get.h_aP(obj)
            %H_AP Addendum, [mm]
            h_aP = 1.0*obj.m;
        end
        
        function h_fP = get.h_fP(obj)
            %H_FP Dedendum, [mm]
            switch obj.type
                case "A"
                    h_fP = 1.25*obj.m;
                case "B"
                    h_fP = 1.25*obj.m;
                case "C"
                    h_fP = 1.25*obj.m;
                case "D"
                    h_fP = 1.40*obj.m;
                otherwise
                    error("prog:input", "Type [%s] is NOT valid.", obj.type);
            end
        end
        
        function h_FfP = get.h_FfP(obj)
            %H_FFP Straight portion of the dedendum, [mm]
            h_FfP = obj.h_fP - obj.c_P;
        end
        
        function h_P = get.h_P(obj)
            %H_P Tooth depth, [mm]
            h_P = obj.h_aP + obj.h_fP;
        end
        
        function h_wP = get.h_wP(obj)
            %H_WP Common depth of rack and tooth, [mm]
            h_wP = obj.h_P - obj.c_P;
        end
        
        function p = get.p(obj)
            %P Pitch, [mm]
            p = pi*obj.m;
        end
        
        function s_P = get.s_P(obj)
            %S_P Tooth thickness, [mm]
            s_P = obj.e_P;
        end
        
        function rho_fP = get.rho_fP(obj)
            %RHO_FP Fillet radius, [mm]
            switch obj.type
                case "A"
                    rho_fP = 0.38*obj.m;
                case "B"
                    rho_fP = 0.30*obj.m;
                case "C"
                    rho_fP = 0.25*obj.m;
                case "D"
                    rho_fP = 0.39*obj.m;
                otherwise
                    error("prog:input", "Type [%s] is NOT valid.", obj.type);
            end
        end
        
    end
end

function mustBeInRange(a, b)
    mustBeGreaterThanOrEqual(a(:), min(b));
    mustBeLessThanOrEqual(   a(:), max(b));
    
%     if any(a(:) < b(1)) || any(a(:) > b(2))
%         error("prog:input", "Value assigned to property is not in range [%f %f]", num2str(b(1)) ,num2str(b(2)));
%     end
end
