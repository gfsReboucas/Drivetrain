classdef Material
    %MATERIAL
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
    
    properties(Constant)
        E         (1, 1) = 206.0e9;  % [Pa],     Young's modulus
        nu        (1, 1) = 0.3;      % [-],      Poisson's ratio
        sigma_Hlim(1, 1) = 1500.0e6; % [Pa],     Allowable contact stress number
        rho       (1, 1) = 7.83e3;   % [kg/m^3], Density
        S_ut      (1, 1) = 700.0e6;  % [Pa],     Tensile strength
        S_y       (1, 1) = 490.0e6;  % [Pa],     Yield strength
    end
    
    methods(Static)
        function disp()
            tab_set = {"Young's modulus",                 "E",            Material.E,          "Pa";
                       "Poisson's ratio",                 "nu",           Material.nu,         "-";
                       "Density",                         "rho",          Material.rho,        "kg/m^3";
                       "Allowable contact stress number", "sigma_Hlim",   Material.sigma_Hlim, "Pa";
                       "Tensile strength",                "S_ut",         Material.S_ut,       "Pa";
                       "Yield strength",                  "S_y",          Material.S_y,        "Pa";
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
    end
end
