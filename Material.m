classdef Material
    %MATERIAL
    % References:
    %   [1] ISO 6336-1:2006 Calculation of load capacity of spur and 
    % helical gears -- Part 1: Basic principles, introduction and general 
    % influence factors.
    %   [2] ISO 6336-2:2006 Calculation of load capacity of spur and helical
    % gears -- Part 2: Calculation of surface durability (pitting)
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
        row = 1; % [-], cf. ISO 6336-2, Table 2
        label string; % 
        E;            % [N/mm^2],  Young's modulus
        nu;           % [-],       Poisson's ratio
        sigma_Hlim;   % [kg/mm^3], Allowable stress number (contact)
        rho;          % [kg/mm^3], Density
        S_ut          % [Pa],      Tensile strength
        S_y;          % [Pa],      Yield strength
        KISS_ID;      % [-],       ID number at KISSsoft's database
    end
    
    properties(Dependent)
        G; % [kg/mm^3], Shear modulus
    end
    
    methods
        function obj = Material(varargin)
            default = {'row'       , 2, ...
                       'label'     , '16MnCr5', ...
                       'E'         , 206.0e3, ...
                       'nu'        , 0.3, ...
                       'sigma_Hlim', 1500.0, ...
                       'rho'       , 7.83e-6, ...
                       'S_ut'      , 700.0e6, ...
                       'S_y'       , 490.0e6, ...
                       'KISS_ID'   , 10250};
                   
            default = scaling_factor.process_varargin(default, varargin);
            
            obj.row        = default.row;
            obj.label      = default.label;
            obj.E          = default.E;
            obj.nu         = default.nu;
            obj.sigma_Hlim = default.sigma_Hlim;
            obj.rho        = default.rho;
            obj.S_ut       = default.S_ut;
            obj.S_y        = default.S_y;
            obj.KISS_ID    = default.KISS_ID;
            
        end
    
        function disp(obj)
            tab_set = {"Label",                             "-",          obj.label,      "-";
                       "Young's modulus",                   "E",          obj.E,          "N/mm^2";
                       "Poisson's ratio",                   "nu",         obj.nu,         "-";
                       "Density",                           "rho",        obj.rho,        "kg/mm^3";
                       "Allowable stress number (contact)", "sigma_Hlim", obj.sigma_Hlim, "N/mm^2";
                       "Tensile strength",                  "S_ut",       obj.S_ut,       "Pa";
                       "Yield strength",                    "S_y",        obj.S_y,        "Pa";
                        };

            Property = tab_set(:, 1);
            Symbol   = tab_set(:, 2);
            Value    = tab_set(:, 3);
            Unit     = tab_set(:, 4);
            
            tab = table(Property, Symbol, Value, Unit, ...
                        'variableNames', ["Property", "Symbol", "Value", "Unit"]);
            
            if(nargout == 0)
                disp(tab);
                clear tab;
            end
        end
        
        function data = export2struct(obj)
            warning('off', 'MATLAB:structOnObject');
            data = struct(obj);
            warning('on', 'MATLAB:structOnObject');
        end
        
    end
    
    %% Get methods:
    methods
        function val = get.G(obj)
            val = (obj.E/2.0)/(1.0 + obj.nu);
        end
    end
end
