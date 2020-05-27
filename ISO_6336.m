classdef ISO_6336
    %ISO_6336 do some calculations following the ISO 6336 standard.
    %
    % References:
    %   [1] ISO 6336-1:2006 Calculation of load capacity of spur and 
    % helical gears -- Part 1: Basic principles, introduction and general 
    % influence factors.
    %   [2] ISO 6336-2:2006 Calculation of load capacity of spur and 
    % helical gears -- Part 2: Calculation of surface durability (pitting).
    %   [3] ISO/TR 6336-30:2017 Calculation of load capacity of spur and 
    % helical gears -- Calculation examples for the application of ISO 6336 
    % parts 1, 2, 3, 5
    %   [4] Arnaudov, K., Karaivanov, D. (2019). Planetary Gear Trains. 
    % Boca Raton: CRC Press, https://doi.org/10.1201/9780429458521
    %
    
    properties(Access = public)
        P;                        % [kW], Rated power
        n_out;                    % [1/min.], Output speed
        gear_set (1, :) Gear_Set; % [-], Gear stage
        S_Hmin;                   % [-], Minimum required safety factor for surface durability (pitting)
        S_Fmin;                   % [-], Minimum required safety factor for tooth bending strength
        L_h;                      % [h], Required life
        K_A;                      % [-], Application factor
    end
    
    properties(Dependent)
        S_H;     % [-], Pitting safety factor
        S_F;     % [-], Tooth bending safety factor
        K_gamma; % [-], Mesh load factor
    end
    
    methods
        function obj = ISO_6336(gset, SHm, SFm, Lh, KA)
            obj.gear_set = gset;
            obj.S_Hmin   = SHm;
            obj.S_Fmin   = SFm;
            obj.L_h      = Lh;
            obj.K_A      = KA;
            
        end
    end
    %% Calculation methods:
    methods
        function SH = Pitting(obj, R_a, R_z, nu_40, line, C_a)
            % preparatory calculations:
            T_1 = (obj.P*1.0e3)/(obj.n_out*pi/30.0);
            
            T_1 = abs(T_1);
            
            SH = obj.calculate_SH(T_1, C_a, E, nu, rho, 
                                  sigma_Hlim, line, nu_40, R_z);
        end
    end
    
    %% Get methods:
    methods
        function val = get.K_gamma(obj)
            Np = obj.gear_set.N_p;
            if(Np == 3)
                val = 1.1;
            elseif(Np == 4)
                val = 1.25;
            elseif(Np == 5)
                val = 1.35;
            elseif(Np == 6)
                val = 1.44;
            elseif(Np == 7)
                val = 1.47;
            else
                val = 1.0;
            end
        end
    end
end
