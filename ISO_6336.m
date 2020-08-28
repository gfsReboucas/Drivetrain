classdef (Abstract) ISO_6336 < Gear_Set
    %ISO_6336 do some calculations following the ISO 6336 standard.
    % see also: KSISO_6336, MATISO_6336, Gear_Set
    %
    % References:
    %   [1] ISO 6336-1:2006 Calculation of load capacity of spur and 
    % helical gears -- Part 1: Basic principles, introduction and general 
    % influence factors.
    %   [2] ISO 6336-2:2006 Calculation of load capacity of spur and 
    % helical gears -- Part 2: Calculation of surface durability (pitting).
    %   [3] ISO 6336-6:2019 Calculation of load capacity of spur and 
    % helical gears -- Part 6: Calculation of service life under variable 
    % load
    %   [4] ISO/TR 6336-30:2017 Calculation of load capacity of spur and 
    % helical gears -- Calculation examples for the application of ISO 6336 
    % parts 1, 2, 3, 5
    %   [5] Arnaudov, K., Karaivanov, D. (2019). Planetary Gear Trains. 
    % Boca Raton: CRC Press, https://doi.org/10.1201/9780429458521
    %
    
    properties(Access = public)
        P_rated;     % [kW],     Rated power
        S_Hmin;      % [-],      Minimum required safety factor for surface durability (pitting)
        S_Fmin;      % [-],      Minimum required safety factor for tooth bending strength
        L_h;         % [h],      Required life
        K_A;         % [-],      Application factor
        n_nominal;   % [1/min.], Nominal speed of each component(gears and carrier)
    end
    
    properties(Access = private)
        idx_fixed; % [-],        Index of the fixed element (if planetary)
        idx_input; % [-],        Index of the input element
        gear_set;
    end
    
    properties(Dependent, Access = public)
        speed_ratio;   % [-],      Speed ratio w.r.t. input element
        torque_ratio;  % [-],      Torque ratio w.r.t. input element
        T_nominal;     % [N-m],    Nominal torque applied on each element
        T_input;       % [N-m],    Input torque of the Gear_Set
        T_output;      % [N-m],    Output torque of the Gear_Set
        n_input;       % [1/min.], Input speed of the Gear_Set
        n_output;      % [1/min.], Output speed of the Gear_Set
        sigma_HP_ref;  % [N/mm^2], Permissive contact stress (reference)
        sigma_HP_stat; % [N/mm^2], Permissive contact stress (static)
        K_gamma;       % [-],      Mesh load factor
        F_tH;          % [N],      Determinant tangential load in a transverse plane for K_Halpha and K_Falpha
    end
    
    properties(Dependent, Access = private)
        idx_output; % [-],        Index of the output element
    end
    
    methods
        function obj = ISO_6336(gset, varargin)
            if(~exist('gset', "var"))
                gset = Gear_Set();
            elseif(~strcmp(class(gset), "Gear_Set"))
                error('gset is not a Gear_Set object.');
            end
            
            T_1 = 9.0e3;
            n_1 = 360.0;
            P_r = 1.0e-3*T_1*n_1*pi/30.0;
            
            default = {'P_rated'    , P_r, ...
                       'S_Hmin'     , 1.0, ...
                       'S_Fmin'     , 1.0, ...
                       'L_h'        , 50.0e3, ...
                       'K_A'        , 1.0, ...
                       'n_nominal'  , [n_1, NaN]};
            
            default = process_varargin(default, varargin);
            
            obj@Gear_Set('configuration', gset.configuration, ...
                         'm_n'          , gset.m_n, ...
                         'alpha_n'      , gset.alpha_n, ...
                         'z'            , gset.z, ...
                         'b'            , gset.b, ...
                         'x'            , gset.x, ...
                         'beta'         , gset.beta, ...
                         'k'            , gset.k, ...
                         'bore_ratio'   , gset.bore_ratio, ...
                         'N_p'          , gset.N_p, ...
                         'a_w'          , gset.a_w, ...
                         'rack_type'    , gset.type, ...
                         'bearing'      , gset.bearing, ...
                         'shaft'        , gset.output_shaft, ...
                         'Q'            , gset.Q, ...
                         'R_a'          , gset.R_a, ...
                         'material'     , gset.material);
            
            obj.gear_set = gset;
            
            obj.P_rated = default.P_rated;
            obj.S_Hmin  = default.S_Hmin;
            obj.S_Fmin  = default.S_Fmin;
            obj.L_h     = default.L_h;
            obj.K_A     = default.K_A;

            [obj.n_nominal, ...
                obj.idx_fixed, obj.idx_input] = obj.set_gear_speed(default.n_nominal);
        end
        
        function [N, sigma] = Pitting_SN_curve(obj)
            N = logspace(4, 10, 100);
            ZN = zeros(2, length(N));
            
            for idx = 1:length(N)
                ZN(:, idx) = obj.life_factor(N(idx));
            end
            
            sigma = diag(obj.sigma_HP_ref)*ZN;
        end
        
        function disp(obj)
            obj.Z_NT
            % to be done...
        end
        
    end
    
    %% Calculation methods:
    methods
%         function variable_load(obj, input_torque, fs)
%             torque = obj.torque_ratio'.*input_torque;
%             
%             torque = torque(:, obj.idx_input);
%             
%             [~, edge_torque, ~] = histcounts(torque);
%             diff_torque = mean(diff(edge_torque));
%             torque_0   = edge_torque(1)   - diff_torque;
%             torque_end = edge_torque(end) + diff_torque;
%             edge_torque = [torque_0, edge_torque, torque_end];
%             
%             [bin_torque, edge_torque, ~] = histcounts(torque, edge_torque);
%             bin_torque = fliplr(bin_torque);
%             edge_torque = fliplr(edge_torque);
%             bin_number = length(bin_torque);
%             
%             edge_speed = obj.P_rated*1.0e3./edge_torque;
%             edge_speed = edge_speed*(30.0/pi);
%             
%             if(strcmp(obj.configuration, 'parallel'))
%                 n_nom = NaN(2, 1);
%             elseif(strcmp(obj.configuration, 'planetary'))
%                 n_nom = NaN(4, 1);
%                 n_nom(obj.idx_fixed) = 0.0;
%             end
%             
%             id = ["ISO_6336:KV", ...
%                   "ISO_6336:SF", ...
%                   "ISO_6336:KS", ...
%                   "MATLAB:COM:InvalidProgid", ...
%                   "CVX:Renamed"];
% 
%             for idx = 1:length(id)
%                 warning("off", id(idx));
%             end
%             
%             gset = obj.gear_set;
%             
%             sig_H = zeros(2, bin_number + 1);
%             SH = zeros(2, bin_number + 1);
%             ZNT = zeros(2, bin_number + 1);
%             N_fail = zeros(2, bin_number + 1);
%             for idx = 1:(bin_number + 1)
%                 n_nom(obj.idx_input) = edge_speed(idx);
%                 copy_obj = ISO_6336(gset, 'P_rated'    , obj.P_rated    , ...
%                                           'S_Hmin'     , obj.S_Hmin     , ...
%                                           'S_Fmin'     , obj.S_Fmin     , ...
%                                           'L_h'        , obj.L_h        , ...
%                                           'K_A'        , obj.K_A        , ...
%                                           'nu_40'      , obj.nu_40      , ...
%                                           'n_nominal'  , n_nom);
% 
%                 sig_H(:, idx) = copy_obj.sigma_H;
%                 SH(:, idx)    = copy_obj.S_H;
%                 ZNT(:, idx)   = copy_obj.Z_NT;
%                 N_fail(:, idx) = copy_obj.num_cycles_failure();
%             end
%             
% %             sig_H = diag(obj.S_H)*sig_H;
%             
%             n = 60.0*obj.L_h*edge_speed*obj.u;
%             
%             for idx = 1:length(id)
%                 warning('on', id(idx));
%             end
%             
%             [N_x, sig_y] = obj.Pitting_SN_curve();
%             
%             for idx = 1:2
%                 aa = find(sig_y(idx, :) <= min(sig_H(idx, :)), 1, 'last')
%                 bb = find(sig_y(idx, :) >= max(sig_H(idx, :)), 1, 'first')
%                 range = aa:bb
%                 N_x = N_x(range);
%                 sig_y = sig_y(:, range);
%             end
%             
%             figure;
%             for idx = 1:2
%                 subplot(1,2,idx)
%                 semilogx(N_x, sig_y(idx, :), 'r');
%                 hold on;
%                 semilogx(n, sigS(idx, :), 'b');
%             end
%             
%             n1 = interp1(sig_y(1, :), log(N_x), log(sig_H(1, :)));
%             n2 = interp1(sig_y(2, :), log(N_x), log(sig_H(2, :)));
%             
%             one2N = 1:bin_number;
%             
%             tab = table(one2N', edge_torque(1:end-1)', edge_torque(2:end)', bin_torque', ...
%                 'variableNames', ["Bin_no", "min_Torque", "Max_Torque", "Load_Cycles"]);
%         end
%         
        function N_star = num_cycles_failure(obj)
            N_star = zeros(1, 2);
            sig_HPR = obj.sigma_HP_ref;
            
            for idx = 1:2
                func = @(x)(obj.sigma_H(idx).*obj.S_H(idx) - obj.life_factor(x).*sig_HPR(idx));
                x0 = obj.N_L(idx)*2.0;
                [N_star(idx),~] = fsolve(func, x0);
            end
        end
        
        function sigH = sigma_H_spectra(obj, torque)
            equiv_force = abs(2.0e3*(torque./obj.d(1))/obj.N_p)*(obj.u + 1.0)/(obj.d(1)*obj.b*obj.u);
            
            sigH = obj.Z_H   * ...
                   obj.Z_E   * ...
                   obj.Z_eps * ...
                   obj.Z_beta* ...
                   obj.Z_BD  * ...
                   sqrt(equiv_force*obj.K_gamma* ...
                                    obj.K_v    * ...
                                    obj.K_Hbeta* ...
                                    obj.K_Halpha);
        end
    end
    
    methods(Static)
        function calc = example_01()
            %EXAMPLE_01 runs the first calculation example from 
            %                      Parameter                         Symbol         Problem       Reference         Calculated         Rel_Error_pc
            %     ___________________________________________    _______________    _______    _______________    _______________    ________________
            %     'Minimum safety factor (pitting)'              'S_Hmin'            ""        [  1.0000e+000]    [  1.0000e+000]    [   0.0000e+000]
            %     'Pitting safety factor for pinion'             'S_H1'              "YES"     [  1.0285e+000]    [  1.0157e+000]    [   1.2505e+000]
            %     'Pitting safety factor for wheel'              'S_H2'              "YES"     [  1.0870e+000]    [  1.0724e+000]    [   1.3390e+000]
            %     'Permissible contact stress for pinion'        'sigma_HP1'         ""        [  1.3385e+003]    [  1.3428e+003]    [-324.6952e-003]
            %     'Permissible contact stress for wheel'         'sigma_HP2'         ""        [  1.4145e+003]    [  1.4178e+003]    [-234.2749e-003]
            %     'Nominal contact stress'                       'sigma_H0'          ""        [  1.2066e+003]    [  1.2065e+003]    [   8.1404e-003]
            %     'Contact stress for pinion'                    'sigma_H1'          "YES"     [  1.3014e+003]    [  1.3221e+003]    [  -1.5951e+000]
            %     'Contact stress for wheel'                     'sigma_H2'          "YES"     [  1.3014e+003]    [  1.3221e+003]    [  -1.5951e+000]
            %     'Nominal tangential load, [N]'                 'F_t'               ""        [127.3520e+003]    [127.3524e+003]    [-299.5607e-006]
            %     'Pitch line velocity, [m/s]'                   'v'                 ""        [  2.6640e+000]    [  2.6642e+000]    [  -7.4461e-003]
            %     'Dynamic factor'                               'K_v'               ""        [  1.0030e+000]    [  1.0005e+000]    [ 246.6092e-003]
            %     'Transverse load factor (root stress)'         'K_Falpha'          ""        [  1.0000e+000]    [  1.0000e+000]    [   0.0000e+000]
            %     'Transverse load factor (contact stress)'      'K_Halpha'          ""        [  1.0000e+000]    [  1.0000e+000]    [   0.0000e+000]
            %     'Face load factor (root stress)'               'K_Fbeta'           "YES"     [  1.1280e+000]    [  1.1601e+000]    [  -2.8446e+000]
            %     'Face load factor (contact stress)'            'K_Hbeta'           "YES"     [  1.1600e+000]    [  1.2002e+000]    [  -3.4681e+000]
            %     'Virtual number of teeth for pinion'           'z_n1'              ""        [ 18.9050e+000]    [ 18.9051e+000]    [-650.6211e-006]
            %     'Virtual number of teeth for wheel'            'z_n2'              ""        [114.5430e+000]    [114.5428e+000]    [ 171.0638e-006]
            %     'Number of load cycles for pinion'             'N_L1'              ""        [  1.0800e+009]    [  1.0800e+009]    [   0.0000e+000]
            %     'Number of load cycles for wheel'              'N_L2'              ""        [178.3000e+006]    [178.2524e+006]    [  26.6813e-003]
            %     'Life factor for pinion'                       'Z_NT1'             ""        [910.0000e-003]    [913.0094e-003]    [-330.6993e-003]
            %     'Life factor for wheel'                        'Z_NT2'             ""        [962.0000e-003]    [964.0118e-003]    [-209.1318e-003]
            %     'Size factor'                                  'Z_X'               ""        [  1.0000e+000]    [  1.0000e+000]    [   0.0000e+000]
            %     'Work hardening factor'                        'Z_W'               ""        [  1.0000e+000]    [  1.0000e+000]    [   0.0000e+000]
            %     'Roughness factor'                             'Z_R'               ""        [965.9900e-003]    [965.9878e-003]    [ 232.3318e-006]
            %     'Velocity factor'                              'Z_v'               ""        [969.1100e-003]    [969.1142e-003]    [-434.0740e-006]
            %     'Lubricant factor'                             'Z_L'               ""        [  1.0474e+000]    [  1.0474e+000]    [ 368.0474e-006]
            %     'Helix angle factor'                           'Z_beta'            ""        [  1.0194e+000]    [  1.0194e+000]    [-366.7828e-006]
            %     'Contact ratio factor'                         'Z_eps'             ""        [803.0000e-003]    [803.3898e-003]    [ -48.5416e-003]
            %     'Elasticity factor'                            'Z_E'               ""        [189.8117e+000]    [189.8117e+000]    [-230.5266e-009]
            %     'Single pair tooth contact factor'             'Z_B'               ""        [  1.0000e+000]    [  1.0000e+000]    [   0.0000e+000]
            %     'Single pair tooth contact factor'             'Z_D'               ""        [  1.0000e+000]    [  1.0000e+000]    [   0.0000e+000]
            %     'Zone factor'                                  'Z_H'               ""        [  2.3953e+000]    [  2.3953e+000]    [-179.4571e-006]
            %     'Single stiffness, [N/(mm-um)]'                'c''                ""        [ 12.3705e+000]    [ 12.3620e+000]    [  68.2642e-003]
            %     'Theoretical single stiffness, [N/(mm-um)]'    'c_th''             ""        [ 17.8558e+000]    [ 17.8436e+000]    [  68.2770e-003]
            %     'Mesh stiffness, [N/(mm-um)]'                  'c_gamma_alpha'     ""        [ 17.4648e+000]    [ 17.4553e+000]    [  54.8918e-003]
            %     'Mesh stiffness, [N/(mm-um)]'                  'c_gamma_beta'      ""        [ 14.8451e+000]    [ 14.8370e+000]    [  54.8750e-003]
            %
            
            mat = Material('row', 2);
            mat = repmat(mat, 1, 2);
            
            gset = Gear_Set('configuration', 'parallel', ...   % configuration
                            'm_n'          , 8.0, ...          % normal module
                            'alpha_n'      , 20.0, ...         % pressure angle
                            'z'            , [17 103], ...     % number of teeth
                            'b'            , 100.0, ...        % face width
                            'x'            , [0.145 0.0], ...  % profile shift coefficient
                            'beta'         , 15.8, ...         % helix angle
                            'k'            , [1 1]*0, ...      % k
                            'bore_ratio'   , [1 1]*0.5, ...    % bore ratio
                            'N_p'          , 1, ...            % number of planets
                            'a_w'          , 500.0, ...        % center distance
                            'rack_type'    , 'D', ...          % rack type
                            'bearing'      , Bearing(), ...    %
                            'shaft'        , Shaft(), ...      %
                            'Q'            , 5.0, ...          % ISO accuracy grade
                            'R_a'          , 1.0, ...          % surface roughness flank
                            'material'     , mat);
            
            T_1 = 9.0e3;
            n_1 = 360.0;
            P_r = 1.0e-3*T_1*n_1*pi/30.0;
            
            calc = ISO_6336(gset, 'P_rated'    , P_r, ...
                                  'S_Hmin'     , 1.0, ...
                                  'S_Fmin'     , 1.0, ...
                                  'L_h'        , 50.0e3, ...
                                  'K_A'        , 1.0, ...
                                  'n_nominal'  , [n_1 NaN], ...
                                  'nu_40'      , 320.0, ...
                                  'C_a'        , 70.0);

            tab_set = {'Minimum safety factor (pitting)'          , 'S_Hmin'       , calc.S_Hmin       , 1.0;
                       'Pitting safety factor for pinion'         , 'S_H1'         , calc.S_H(1)       , 1.02853;
                       'Pitting safety factor for wheel'          , 'S_H2'         , calc.S_H(2)       , 1.08696;
                       'Permissible contact stress for pinion'    , 'sigma_HP1'    , calc.sigma_HP(1)  , 1338.48050;
                       'Permissible contact stress for wheel'     , 'sigma_HP2'    , calc.sigma_HP(2)  , 1414.52551;
                       'Nominal contact stress'                   , 'sigma_H0'     , calc.sigma_H0     , 1206.58207;
                       'Contact stress for pinion'                , 'sigma_H1'     , calc.sigma_H(1)   , 1301.35343;
                       'Contact stress for wheel'                 , 'sigma_H2'     , calc.sigma_H(2)   , 1301.35343;
                       'Nominal tangential load, [N]'             , 'F_t'          , calc.F_t          , 127352.0;
                       'Pitch line velocity, [m/s]'               , 'v'            , calc.v_pitch_line , 2.664;
                       'Dynamic factor'                           , 'K_v'          , calc.K_v          , 1.003;
                       'Transverse load factor (root stress)'     , 'K_Falpha'     , calc.K_Falpha     , 1.0;
                       'Transverse load factor (contact stress)'  , 'K_Halpha'     , calc.K_Halpha     , 1.0;
                       'Face load factor (root stress)'           , 'K_Fbeta'      , calc.K_Fbeta      , 1.12803;
                       'Face load factor (contact stress)'        , 'K_Hbeta'      , calc.K_Hbeta      , 1.16;
                       'Virtual number of teeth for pinion'       , 'z_n1'         , calc.z_n(1)       , 18.905;
                       'Virtual number of teeth for wheel'        , 'z_n2'         , calc.z_n(2)       , 114.543;
                       'Number of load cycles for pinion'         , 'N_L1'         , calc.N_L(1)       , 1.080e9;
                       'Number of load cycles for wheel'          , 'N_L2'         , calc.N_L(2)       , 1.783e8;
                       'Life factor for pinion'                   , 'Z_NT1'        , calc.Z_NT(1)      , 0.91;
                       'Life factor for wheel'                    , 'Z_NT2'        , calc.Z_NT(2)      , 0.962;
                       'Size factor'                              , 'Z_X'          , calc.Z_X          , 1.0;
                       'Work hardening factor'                    , 'Z_W'          , calc.Z_W          , 1.0;
                       'Roughness factor'                         , 'Z_R'          , calc.Z_R          , 0.96599;
                       'Velocity factor'                          , 'Z_v'          , calc.Z_v          , 0.96911;
                       'Lubricant factor'                         , 'Z_L'          , calc.Z_L          , 1.04739;
                       'Helix angle factor'                       , 'Z_beta'       , calc.Z_beta       , 1.01944;
                       'Contact ratio factor'                     , 'Z_eps'        , calc.Z_eps        , 0.803;
                       'Elasticity factor'                        , 'Z_E'          , calc.Z_E          , 189.81170;
                       'Single pair tooth contact factor'         , 'Z_B'          , calc.Z_BD(1)      , 1.0;
                       'Single pair tooth contact factor'         , 'Z_D'          , calc.Z_BD(2)      , 1.0;
                       'Zone factor'                              , 'Z_H'          , calc.Z_H          , 2.39533;
                       'Single stiffness, [N/(mm-um)]'            , 'c'''          , calc.cprime       , 12.37047;
                       'Theoretical single stiffness, [N/(mm-um)]', 'c_th'''       , calc.cprime_th    , 17.85584;
                       'Mesh stiffness, [N/(mm-um)]'              , 'c_gamma_alpha', calc.c_gamma_alpha, 17.46485;
                       'Mesh stiffness, [N/(mm-um)]'              , 'c_gamma_beta' , calc.c_gamma_beta , 14.84512;
                       };

            Parameter = tab_set(:, 1);
            Symbol    = tab_set(:, 2);
            val_calc  = tab_set(:, 3);
            val_ref   = tab_set(:, 4);
            rel_err   = (1.0 - cell2mat(val_calc)./cell2mat(val_ref))*100.0;
            big_err = repmat("", length(rel_err), 1);
            big_err(abs(rel_err) > 1.0) = 'YES';
            rel_err = mat2cell(rel_err, ones(length(rel_err),1));
            
            tab = table(Parameter, Symbol, big_err, val_ref, val_calc, rel_err, ...
                'variableNames', ["Parameter", "Symbol", "Problem", "Reference", "Calculated", "Rel_Error_pc"]);
            disp(tab)
        end
        
        function calc = example_02()
            mat = Material('row', 2);
            mat = repmat(mat, 1, 2);

            gset = Gear_Set('configuration', 'parallel'             , ... % configuration
                            'm_n'          , 8.467                  , ... % normal module
                            'alpha_n'      , 25.0                   , ... % pressure angle
                            'z'            , [17 60]                , ... % number of teeth
                            'b'            , 152.4                  , ... % face width
                            'x'            , [0.172 0.0015]         , ... % profile shift coefficient
                            'beta'         , 15.5                   , ... % helix angle
                            'k'            , [1 1]*0                , ... % k
                            'bore_ratio'   , [1 1]*0.5              , ... % bore ratio
                            'N_p'          , 1                      , ... % number of planets
                            'a_w'          , 339.727                , ... % center distance
                            'rack_type'    , 'D'                    , ... % rack type
                            'bearing'      , repmat(Bearing(), 1, 6), ... %
                            'shaft'        , Shaft()                , ... %
                            'Q'            , 6.0                    , ... % ISO accuracy grade
                            'R_a'          , 1.0                    , ... % surface roughness flank
                            'material'     , mat);

            T_1 = mean([25347 25423]); % [N-m], from [6], Table 4, bin 3
            n_1 = 35.2;                % [1.0/min.]
            P_r = 1.0e-3*T_1*n_1*pi/30.0;
            
            calc = ISO_6336(gset, 'P_rated'    , P_r      , ...
                                  'S_Hmin'     , 1.0      , ...
                                  'S_Fmin'     , 1.0         , ...
                                  'L_h'        , 1.0e3, ...
                                  'K_A'        , 1.0      , ...
                                  'n_nominal'  , [n_1 NaN], ...
                                  'nu_40'      , 320.0    , ...
                                  'C_a'        , 70.0);

        end
        
        function [N, edges, bin] = LDD(signal, varargin)
            % adapted from Simpack's Load Calculation Tool [1], script 
            % SC_ldd.qs
            %
            % [1] Simpack Wind. Access on 28.08.2020. 
            % https://www.3ds.com/products-services/simulia/training/course-descriptions/simpack-wind/
            %
            
            num_bins = 40;
            LDD_data = zeros(1, num_bins);
            for idx = 1:length(signal)
                delta = signal(idx) - bin_min;
                act_test = floor(delta/bin_step);
                
                if(act_test == -1)
                    bin_test = 1;
%                 elseif(bin_test == num_bins)
%                     bin_test = num_bins - 1;
                end
                % oc is the occurrence?
                LDD_data(bin_test) = LDD_data(bin_test) + time_step*oc;
            end
            
%             [x_min, x_Max] = bounds(signal);
            
%             edges = linspace(x_min, x_Max, 200);
%             default = {};
            [N,edges,bin] = histcounts(signal);
        end
    end

    methods(Access = private)
        function val = F_t_spectra(obj, torque)
            val = 2.0e3*(torque/obj.d(1))/obj.N_p;
            val = abs(val);
        end
        
        function val = F_tH_spectra(obj, torque)
            val = 2.0e3*(torque/obj.d(1))/obj.N_p;
            val = abs(val);
        end
        
        function [n_new, idx_zero, idx_inp] = set_gear_speed(obj, n_old)
            %SET_GEAR_SPEED sets the speeds of a Gear_Set. The array n_old
            % should specify only the fixed element (if a planetary
            % Gear_Set) and the input element. The other speeds are
            % calculated in this method and should be defined as NaN in
            % n_old.
            %
            
            if(strcmpi(obj.configuration, 'parallel'))
                idx_zero = nan;
                
                U = diag([1.0/obj.u, obj.u]);
                U = flip(U);
            elseif(strcmpi(obj.configuration, 'planetary'))
                idx_zero = find(n_old == 0);
                n_old(idx_zero) = nan;
                zc = num2cell(obj.z);
                [s, p, r] = deal(zc{:});
                r = abs(r);
                U = zeros(4);

                switch(idx_zero)
                    case 1 % fixed: sun
                           % input: ring*
                           % output: carrier* (* or vice-versa)
                        U(:, 4) = [0.0, (r + s)/p, (1.0 + s/r), 1.0]';
                    case 2
                        error('ISO_6336:speed', 'Direct drive not implemented yet.');
                    case 3 % fixed: ring
                           % input: sun*
                           % output: carrier*
                        U(:, 4) = [(1.0 + r/s), (r + s)/p, 0.0, 1.0]';
                    case 4 % fixed: carrier
                           % input: sun*
                           % output: ring*
                        U(:, 3) = [-r/s, r/p, 1.0, 0.0];
                    otherwise
                        if(isempty(idx_zero))
                            error('ISO_6336:speed', 'There are no fixed elements or the size.');
                        elseif(idx_zero > 4)
                            error('ISO_6336:speed', 'There should be a maximum of 4 velocities.');
                        end
                end
                
            end
            
            A = eye(length(U)) - U;
            x = null(A);
            
            idx_inp = find(~isnan(n_old));
            x = x./x(idx_inp);
            n_new = x*n_old(idx_inp);
        end
        
        function val = dynamic_factor(obj)
            z1 = obj.z(1);
            uu = obj.u;
            cond = (obj.v_pitch_line*z1/100.0)*sqrt((uu.^2)/(1.0 + uu.^2));
            if(cond < 3.0) % [m/s]
                warning('ISO_6336:KV', ['Calculating K_v using method B', ...
                    ' outside of its useful range. ', ...
                'More info at the end of Sec. 6.3.2 of ISO 6336-1.']);
            end
            
            rho = [obj.material.rho];
            
            % transverse base pitch deviation, [um]:
            fpb = obj.f_pb;
            
            % profile form deviation, [um]
            f_falpha = max(obj.f_falpha);
            
            % Estimated running allowance (pitch deviation):
            y_p = obj.y_alpha;
            
            % Estimated running allowance (flank deviation):
            y_f = f_falpha*75.0e-3;
            
            % transverse effective base pitch deviation after running-in, [um]:
            f_pbeff = fpb - y_p;
            
            % effective profile form deviation after running-in, [um]
            f_falphaeff = f_falpha - y_f; % 0.925*y_f
            
            if(strcmp(obj.configuration, 'parallel'))
                % Based on Sec. 6.5.9 of ISO 6336-1, Eq. (30), assuming gears:
                % - of solid construction, and
                % - with the same density
                
                num = pi*rho(1)*(obj.u^2)*obj.d_m(1)^4;
                den = 8.0*(obj.u^2 + 1.0)*obj.d_b(1)^2;
                m_red = num/den;
                
                % Resonance running speed, Eq. (6) [1]:
                n_E1 = 3.0e4*sqrt(obj.c_gamma_alpha/m_red)/(pi*obj.z(1)); % [1/min]
                
                % Resonance ratio, Eq. (9) [1]:
                N = obj.n_output/n_E1;
                
                val = obj.dynamic_factor_from_range(N, obj.C_a, ...
                    f_pbeff, f_falphaeff);
            elseif(strcmp(obj.configuration, 'planetary'))
                m_sun = rho(1)*(pi/8.0)*(obj.d_m(1)^4)/(obj.d_b(1)^2); % sun
                m_pla = rho(2)*(pi/8.0)*(obj.d_m(2)^4)/(obj.d_b(2)^2); % planet
                m_red1 = (m_pla*m_sun)/(m_sun + m_pla*obj.N_p);
                m_red2 =  m_pla;
                
                % Resonance running speed:
                n_E11 = 3.0e4*sqrt(obj.cprime/m_red1)*1.0/(pi*obj.z(1)); % [1/min]
                n_E12 = 3.0e4*sqrt(obj.cprime/m_red2)*1.0/(pi*obj.z(2)); % [1/min]
                
                % Resonance ratio:
                N_1 =  obj.n_output/n_E11;
                N_2 = (obj.n_output/obj.u)/n_E12;

                K_v1 = obj.dynamic_factor_from_range(N_1, obj.C_a, ...
                    f_pbeff, f_falphaeff);
                
                K_v2 = obj.dynamic_factor_from_range(N_2, obj.C_a, ...
                    f_pbeff, f_falphaeff);
                
                val = max(K_v1, K_v2);
            end
        end
        
        function K_v = dynamic_factor_from_range(obj, N, C_a, f_pbeff, f_falphaeff)
            % Table 8, [1]:
            C_v1 = 0.32;
            C_v5 = 0.47;
            
            eps_g = obj.eps_gamma;
            if((1.0 < eps_g) && (eps_g <= 2.0))
                C_v2 = 0.34;
                C_v3 = 0.23;
                C_v4 = 0.90;
                C_v6 = 0.47;
            elseif(eps_g > 2.0)
                C_v2 =  0.57 /(eps_g - 0.3);
                C_v3 =  0.096/(eps_g - 1.56);
                C_v4 = (0.57 - 0.05*eps_g)/(eps_g - 1.44);
                C_v6 =  0.12/(eps_g - 1.74);
            end
                
            if((1.0 < eps_g) && (eps_g <= 1.5))
                C_v7 = 0.75;
            elseif((1.5 < eps_g) && (eps_g <= 2.5))
                C_v7 = 0.125*sin(pi*(eps_g - 2.0)) + 0.875;
            elseif(eps_g > 2.5)
                C_v7 = 1.0;
            end
            
            cp = obj.cprime;
            B_p = cp*f_pbeff/obj.line_load;        % Eq. (15)
            B_f = cp*f_falphaeff/obj.line_load;    % Eq. (16)
            B_k = abs(1.0 - cp*C_a/obj.line_load); % Eq. (17)
            
            % Dynamic factor:
            if(obj.line_load < 100) % [N/mm]
                N_S = 0.5 + 0.35*sqrt(obj.line_load/100.0); % Eq. (11), [1]
            else
                N_S = 0.85; % Eq. (12) [1]
            end
            
            if(N <= N_S)
                % Section 6.5.4, Eq. (13), [1]:
                K = C_v1*B_p + C_v2*B_f + C_v3*B_k;
                K_v = N*K + 1.0; % Eq. (14)
            elseif((N_S < N) && (N <= 1.15))
                % Section 6.5.5, Eq. (20), [1]:
                K_v = C_v1*B_p + C_v2*B_f + C_v4*B_k + 1.0;
            elseif((1.15 < N) && (N < 1.5))
                K_vN115 = C_v1*B_p + C_v2*B_f + C_v4*B_k + 1.0;
                K_vN15  = C_v5*B_p + C_v6*B_f + C_v7;
                % Section 6.5.7, Eq. (22), [1]:
                K_v = K_vN15 + (K_vN115 - K_vN15)*(1.5 - N)/0.35;
            elseif(N >= 1.5)
                % Section 6.5.6, Eq. (21), [1]:
                K_v = C_v5*B_p + C_v6*B_f + C_v7;
            end
        end
        
        function [Z_L, Z_v] = lubrication_velocity_factor(obj)

            % [N/mm^2],  Allowable contact stress number:
            sig_Hlim = min(obj.sigma_Hlim);
            
            if(sig_Hlim  < 850.0) % [N/mm^2]
                C_ZL = 0.83;
            elseif((850.0 <= sig_Hlim) && (sig_Hlim  < 1200.0))
                C_ZL = sig_Hlim/4375.0 + 0.6357;
            else
                C_ZL = 0.91;
            end
            
            % Lubricant factor:
            Z_L = C_ZL + 4.0*(1.0 - C_ZL)/(1.2 + 134.0/obj.nu_40)^2;
            
            % Velocity factor:
            C_Zv = C_ZL + 0.02;
            Z_v = C_Zv + 2.0*(1.0 - C_Zv)/sqrt(0.8 + 32.0/obj.v_pitch_line);
        end
        
        function ZN = life_factor(obj, N)
            
            ZN = zeros(1, 2);
            
            for idx = 1:2
                line = obj.material(idx).row;
                ratio = obj.Z_NT(idx); % obj.sigma_HP_stat(idx)/obj.sigma_HP_ref(idx);
                
                switch line
                    case 1
                        % St, V, GGG (perl. bai.), GTS (perl.), Eh, IF (when limited pitting is permitted)
                        if(N < 6.0e5)
                            val = power(3.0e8/6.0e5, ...
                                0.3705*log(ratio));
                        elseif((6.0e5 < N) && (N <= 1.0e7))
                            val = power(3.0e8/N, ...
                                0.3705*log(ratio));
                        elseif((1.0e7 < N) && (N <= 1.0e9))
                            val = power(1.0e9/N, ...
                                0.2791*log(ratio));
                        else
                            val = 1.0;
                        end
                    case 2
                        % St, V, GGG (perl. bai.), GTS (perl.), Eh, IF
                        if(N < 1.0e5)
                            val = power(5.0e7/1.0e5, ...
                                0.3705*log(ratio));
                        elseif((1.0e5 < N) && (N <= 5.0e7))
                            val = power(5.0e7/N, ...
                                0.3705*log(ratio));
                        else
                            val = 1.0;
                        end
                    case 3
                        % GG, GGG (ferr.), NT (nitr.), NV (nitr.)
                        if(N < 1.0e5)
                            val = power(2.0e6/1.0e5, ...
                                0.7686*log(ratio));
                        elseif((1.0e5 < N) && (N <= 2.0e6))
                            val = power(2.0e6/N, ...
                                0.7686*log(ratio));
                        else
                            val = 1.0;
                        end
                    case 4
                        % NV (nitrocar.)
                        if(N < 1.0e5)
                            val = power(2.0e6/1.0e5, ...
                                0.7686*log(ratio));
                        elseif((1.0e5 < N) && (N <= 2.0e6))
                            val = power(2.0e6/N, ...
                                0.7686*log(ratio));
                        else
                            val = 1.0;
                        end
                    otherwise
                        error("ISO_6336:ZN", "Invalid input [%d].\nValid options are 1 to 4.\n", line);
                end
                ZN(idx) = val;
            end
            
        end
        
    end
    
    %% Set methods:
    methods
    end
    
    %% Get methods:
    methods
        function val = get.speed_ratio(obj)
            val = obj.n_nominal./obj.n_input;
        end
        
        function val = get.torque_ratio(obj)
            val = obj.T_nominal./obj.T_input;
        end
        
        function val = get.F_tH(obj)
            val = obj.F_t*obj.K_gamma* ...
                          obj.K_A    * ...
                          obj.K_v    * ...
                          obj.K_Hbeta;
        end
        
        function val = get.sigma_HP_ref(obj)
            ZNT = 1.0;
            val = [obj.material(1:2).sigma_Hlim].*ZNT*obj.Z_L* ...
                                                      obj.Z_v* ...
                                                      obj.Z_R* ...
                                                      obj.Z_W* ...
                                                      obj.Z_X./obj.S_Hmin;
        end
        
        function val = get.sigma_HP_stat(obj)
            val = [obj.material(1:2).sigma_Hlim].*obj.Z_NT* ...
                                                  obj.Z_L* ...
                                                  obj.Z_v* ...
                                                  obj.Z_R* ...
                                                  obj.Z_W* ...
                                                  obj.Z_X./obj.S_Hmin;
        end
        
        function val = get.T_nominal(obj)
            val = (obj.P_rated*1.0e3)./(obj.n_nominal*pi/30.0);
            val(2) = val(2)/obj.N_p;
            val(val == inf) = nan;
        end
        
        function val = get.T_input(obj)
            val = obj.T_nominal(obj.idx_input);
        end
        
        function val = get.T_output(obj)
            val = obj.T_nominal(obj.idx_output);
        end
        
        function val = get.n_input(obj)
            val = obj.n_nominal(obj.idx_input);
        end
        
        function val = get.n_output(obj)
            val = obj.n_nominal(obj.idx_output);
        end
        
        function val = get.idx_output(obj)
            if(strcmpi(obj.configuration, 'parallel'))
                val = setdiff([1 2], obj.idx_input);
            elseif(strcmpi(obj.configuration, 'planetary'))
                switch(obj.idx_fixed)
                    case 1
                        val = setdiff([3 4], obj.idx_input);
                    case 2
                        error('ISO_6336:idx_out', 'Direct drive not implemented yet.');
                    case 3
                        val = setdiff([1 4], obj.idx_input);
                    case 4
                        val = setdiff([1 3], obj.idx_input);
                    otherwise
                        error('ISO_6336:idx_out', 'Invalid index.');
                end
            end
        end
        
        function val = get.K_gamma(obj)
            switch(obj.N_p)
                case 3
                    val = 1.1;
                case 4
                    val = 1.25;
                case 5
                    val = 1.35;
                case 6
                    val = 1.44;
                case 7
                    val = 1.47;
                otherwise
                    val = 1.0;
            end
        end
        
    end
end
