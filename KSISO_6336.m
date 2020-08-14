classdef KSISO_6336 < ISO_6336
    %KSISO_6336 do some calculations following the ISO 6336 standard using
    % KISSsoft COM interface.
    % see also: ISO_6336, MATISO_6336, Gear_Set, KISSsoftCOM
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
    
    properties(Access = private)
        ks KISSsoftCOM;
    end
    
    properties(Access = public)
        save_report;
        show_report;
    end
    
    properties(Dependent, Access = public)
        S_H;           % [-],      Pitting safety factor
        S_F;           % [-],      Tooth bending safety factor
        sigma_H0;      % [N/mm^2], Nominal contact stress
        sigma_H;       % [N/mm^2], Contact stress
        sigma_HP;      % [N/mm^2], Permissible contact stress
        sigma_Hlim;    % [N/mm^2], Allowable stress number
        F_t;           % [N],      Transverse tangential load at reference cylinder (nominal, per mesh)
        line_load;     % [N/mm],   Line load
        v_pitch_line;  % [m/s],    Pitch line velocity
        N_L;           % [-],      Number of load cycles
        K_Falpha;   % [-],        Transverse load factor (root stress) 
        K_Fbeta;    % [-],        Face load factor (root stress) 
        K_Halpha;   % [-],        Transverse load factor (contact stress)
        K_Hbeta;    % [-],        Face load factor (contact stress)
        K_v;        % [-],        Dynamic factor
        Z_beta;     % [-],        Helix angle factor
        Z_eps;      % [-],        Contact ratio factor
        Z_BD;       % [-],        Single pair tooth contact factor
        Z_E;        % [N^0.5/mm], Elasticity factor
        Z_H;        % [-],        Zone factor
        Z_L;        % [-],        Lubrication factor
        Z_N;        % [-],        Life factor for contact stress
        Z_NT;       % [-],        Life factor for contact stress for reference test conditions
        Z_R;        % [-],        Roughness factor affecting surface durability
        Z_v;        % [-],        Velocity factor (circumferential velocity at pitch line)
        Z_W;        % [-],        Work hardening factor
        Z_X;        % [-],        Size factor (pitting)
    end
    
    methods
        function obj = KSISO_6336(gset, varargin)
            if(~KISSsoftCOM.is_installed())
                error('KSISO_6336:not_installed', 'KISSsoft\x00a9 COM interface not found.');
            end
            
            if(~exist('gset', "var"))
                gset = Gear_Set();
            elseif(~strcmp(class(gset), "Gear_Set"))
                error('gset is not a Gear_Set object.');
            end
            
            T_1 = 9.0e3;
            n_1 = 360.0;
            P_r = 1.0e-3*T_1*n_1*pi/30.0;
            
            default = {'P_rated'     , P_r, ...
                       'S_Hmin'      , 1.0, ...
                       'S_Fmin'      , 1.0, ...
                       'L_h'         , 50.0e3, ...
                       'K_A'         , 1.0, ...
                       'lubricant_ID', 11170, ...
                       'n_nominal'   , [n_1, NaN]};
            
            default = process_varargin(default, varargin);
            
            obj@ISO_6336(gset, 'P_rated'    , default.P_rated, ...
                               'S_Hmin'     , default.S_Hmin, ...
                               'S_Fmin'     , default.S_Fmin, ...
                               'L_h'        , default.L_h, ...
                               'K_A'        , default.K_A, ...
                               'n_nominal'  , default.n_nominal);
            
            if(strcmpi(obj.configuration, 'parallel'))
                module = 'Z012';
                std_file = 'CylGearPair 1 (spur gear).Z12';
            elseif(strcmpi(obj.configuration, 'planetary'))
                module = 'Z014';
                std_file = 'PlanetarySet 1 (ISO6336).Z14';
            end
            
            switch(obj.type)
                case 'A'
                    ref_prof_ID = 10010;
                case 'B'
                    ref_prof_ID = 10011;
                case 'C'
                    ref_prof_ID = 10012;
                case 'D'
                    ref_prof_ID = 10013;
                otherwise
                    error('Reference profile [%s] not supported.', obj.type);
            end
            
            ksCOM = KISSsoftCOM(module);
            version = ksCOM.get_version();
            version = strrep(version, '/', '-');

            file_name = sprintf('C:\\Program Files (x86)\\KISSsoft %s\\example\\%s', version, std_file);

            try
                ksCOM.load_file(file_name);
            catch err
                delete(ksCOM);
                error(err.identifier, '%s', err.message);
            end
            
            ksCOM.set_var('RechSt.RechenMethID' , 10030); % calculation method: ISO 6336:2006 Method B
%             ksCOM.set_var('RechSt.GeometrieMeth', true);  % tooth geometry according to ISO 21771:2007
            
            ksCOM.set_var('ZS.AnzahlZwi',         obj.N_p     ); % number of planets
            ksCOM.set_var('ZS.Geo.mn'   ,         obj.m_n     ); % normal module
            ksCOM.set_var('ZP[0].a'     ,         obj.a_w     ); % center distance
            ksCOM.set_var('ZS.Geo.alfn' , deg2rad(obj.alpha_n)); % normal pressure angle
            ksCOM.set_var('ZS.Geo.beta' , deg2rad(obj.beta)   ); % helix angle

            for idx = 1:numel(obj.z)
                ksCOM.set_var(sprintf('ZR[%d].Tool.type'           , idx - 1), 2              ); % reference profile gear
                ksCOM.set_var(sprintf('ZR[%d].z'                   , idx - 1), obj.z(idx)     ); % number of teeth
                ksCOM.set_var(sprintf('ZR[%d].x.nul'               , idx - 1), obj.x(idx)     ); % profile shift coefficient
                ksCOM.set_var(sprintf('ZR[%d].b'                   , idx - 1), obj.b          ); % face width
                ksCOM.set_var(sprintf('ZR[%d].Tool.RefProfile.DBID', idx - 1), ref_prof_ID    ); % Rack profile ID
                ksCOM.set_var(sprintf('ZR[%d].Vqual'               , idx - 1), obj.Q          ); % ISO accuracy grade
                ksCOM.set_var(sprintf('ZR[%d].RAH'                 , idx - 1), obj.R_a        ); % Arithmetic mean roughness
                ksCOM.set_var(sprintf('ZR[%d].RAF'                 , idx - 1), obj.R_z        );
                ksCOM.set_var(sprintf('ZR[%d].mat.DBID', idx - 1)  , obj.material(idx).KISS_ID); % material ID
            end
            
            if(strcmpi(obj.configuration, 'planetary'))
                ksCOM.set_var('ZP[0].Impulsnominal', 0); % carrier as driving?
                
                ksCOM.set_var('ZP[1].a', -obj.a_w); % center distance
                ksCOM.set_var('ZP[2].a', 0.0     );
            end
            
            ksCOM.set_var('ZS.P'          , obj.P_rated); % rated power
            ksCOM.set_var('ZS.Pnominal'   , obj.P_rated); % rated power
            ksCOM.set_var('ZR[0].Tnominal', obj.T_nominal(1)); % torque at gear 1
            ksCOM.set_var('ZR[0].n'       , obj.n_nominal(1)); % speed at gear 1

            % Safeties are NOT depending on size
            ksCOM.set_var('ZS.SSiFlags.SizeDependent', false);
            
            ksCOM.set_var('ZS.SSi.Flanke'      , obj.S_Hmin);
            ksCOM.set_var('ZS.SSi_0to05.Flanke', obj.S_Hmin);
            ksCOM.set_var('ZS.SSi_at1.Flanke'  , obj.S_Hmin);
            ksCOM.set_var('ZS.SSi_2to20.Flanke', obj.S_Hmin);
            ksCOM.set_var('ZS.SSi_fix.Flanke'  , obj.S_Hmin);
            
            ksCOM.set_var('ZS.SSi.Fuss'      , obj.S_Fmin);
            ksCOM.set_var('ZS.SSi_0to05.Fuss', obj.S_Fmin);
            ksCOM.set_var('ZS.SSi_at1.Fuss'  , obj.S_Fmin);
            ksCOM.set_var('ZS.SSi_2to20.Fuss', obj.S_Fmin);
            ksCOM.set_var('ZS.SSi_fix.Fuss'  , obj.S_Fmin);
            
            ksCOM.set_var('ZS.H'   , obj.L_h);
            ksCOM.set_var('ZS.KA'  , obj.K_A);
            ksCOM.set_var('ZS.Kgam', obj.K_gamma);
            ksCOM.set_var('ZS.Oil.SchmierTypID', default.lubricant_ID);
            
            if(~ksCOM.check_calculations())
                delete(ksCOM);
                error('Problem with KISSsoft calculation.');
            end
            
            obj.ks = ksCOM;
            obj.ks.show_UI();
        end
        
        function save_KS(obj, name)
            obj.ks.save_file(name);
        end
    end
    
    %% Calculation methods:
    methods
        function [SH, SF] = safety_factors(obj, varargin)
            SH = obj.S_H;
            SF = obj.S_F;
        end
        
        function [SH, SF] = KISSsoft_calculation(obj, varargin)
            default = {'lubricant_ID', '11220', ...
                       'stage_idx'   ,  0, ...
                       'save_report' ,  false, ...
                       'show_report' ,  false};
            
            default = process_varargin(default, varargin);
            
%             ks = obj.KISSsoft();
            ksCOM.set_var('ZS.P'       , obj.P_rated); % rated power
            ksCOM.set_var('ZS.Pnominal', obj.P_rated); % rated power
            
            ksCOM.set_var('ZS.Oil.SchmierTypID', default.lubricant_ID);
            
            flag_impul = true;
            if(strcmp(obj.configuration, "parallel"))
                flag_impul = true;
                template_code = 12;
            elseif(strcmp(obj.configuration, 'planetary'))
                n_c = obj.n_nominal(4);
                T_c = obj.T_nominal(4);
                n_p = obj.n_nominal(2);
                ksCOM.set_var('ZR[1].n'               , abs(n_p)); % speed at gear 2: planet
                ksCOM.set_var('ZR[2].n'               , 0.0     ); % speed at gear 3: ring
                ksCOM.set_var('ZS.Planet.nStegnominal', n_c     ); % speed at planet carrier
                ksCOM.set_var('ZS.Planet.nSteg'       , n_c     ); % speed at planet carrier
                
                ksCOM.set_var('ZR[1].T'               , abs(n_p)); % speed at gear 2: planet
                ksCOM.set_var('ZS.Planet.TStegnominal', T_c     ); % torque at planet carrier
                ksCOM.set_var('ZS.Planet.nSteg'       , T_c     ); % torque at planet carrier
                
                flag_impul = false;
                template_code = 14;
            end
            
            ksCOM.set_var('ZR[0].Tnominal', obj.T_nominal(1)); % torque at gear 1
            ksCOM.set_var('ZR[0].n'       , obj.n_nominal(1)); % speed at gear 1
            ksCOM.set_var('ZP[0].Impuls'  , flag_impul); % sense of rotation, gear 1

            ksCOM.set_var('ZS.SSi.Flanke'      , obj.S_Hmin);
            ksCOM.set_var('ZS.SSi_0to05.Flanke', obj.S_Hmin);
            ksCOM.set_var('ZS.SSi_at1.Flanke'  , obj.S_Hmin);
            ksCOM.set_var('ZS.SSi_2to20.Flanke', obj.S_Hmin);
            ksCOM.set_var('ZS.SSi_fix.Flanke'  , obj.S_Hmin);
            
            ksCOM.set_var('ZS.SSi.Fuss'      , obj.S_Fmin);
            ksCOM.set_var('ZS.SSi_0to05.Fuss', obj.S_Fmin);
            ksCOM.set_var('ZS.SSi_at1.Fuss'  , obj.S_Fmin);
            ksCOM.set_var('ZS.SSi_2to20.Fuss', obj.S_Fmin);
            ksCOM.set_var('ZS.SSi_fix.Fuss'  , obj.S_Fmin);
            
            ksCOM.set_var('ZS.H'   , obj.L_h);
            ksCOM.set_var('ZS.KA'  , obj.K_A);
            ksCOM.set_var('ZS.Kgam', obj.K_gamma);
            
            if(~ksCOM.calculate())
%                 delete(ks);
            end
            
            kv = ksCOM.get_var('ZS.KVcalc');
            
            flag_kv = false;
            if(kv < 1.05)
                flag_kv = true;
                for idx = 1:2
                    ksCOM.set_var(sprintf('ZP[%d].KV.KV', idx), 1.05);
                end
                
                ksCOM.set_var('Zst.KVFlag', flag_kv);
            end
            
            flag_khb = false;
            for idx = 1:2
                khb = ksCOM.get_var(sprintf('ZP[%d].KHb', idx - 1));
                
                if(khb < 1.15)
                    flag_khb = true;
                    ksCOM.set_var(sprintf('ZP[%d].KHb', idx - 1), 1.15);
                end
            end
            
            if(flag_khb)
                ksCOM.set_var('Zst.KHbVariant', true);
            end
            
            if(~ksCOM.calculate())
%                 delete(ks);
            end
            
            if(flag_kv || flag_khb)
                ksCOM.calculate();
            end
            
            max_Np = obj.N_p + 1;
            
            SH = zeros(max_Np, 1);
            SF = zeros(max_Np, 1);
            
            for idx = 1:length(obj.z)
                SH(idx) = ksCOM.get_var(sprintf('ZPP[%d].Flanke.SH', idx - 1));
                SF(idx) = ksCOM.get_var(sprintf('ZPP[%d].Fuss.SF'  , idx - 1));
            end
            
            SH(SH == 0.0) = [];
            SF(SF == 0.0) = [];
            
            SH(isnan(SH)) = [];
            SF(isnan(SF)) = [];
            
            if(isrow(SH))
                SH = SH';
            end
            
            if(isrow(SF))
                SF = SF';
            end
            
            if(default.save_report)
                if(strcmpi(class(obj), "Drivetrain"))
                    report_name = sprintf("%s\\stage_%02d.rtf" , pwd, default.stage_idx);
                    file_name   = sprintf("%s\\stage_%02d.Z0%d", pwd, default.stage_idx, template_code);
                else
                    report_name = sprintf("%s\\@%s\\stage_%02d.rtf" , pwd, class(obj), default.stage_idx);
                    file_name   = sprintf("%s\\@%s\\stage_%02d.Z0%d", pwd, class(obj), default.stage_idx, template_code);
                end
                
                ksCOM.save_file(file_name);
                version = ksCOM.get_version;
                template = sprintf("C:\\Program Files (x86)\\KISSsoft %s\\rpt\\Z0%dLe0.RPT", version, template_code);
                
                ksCOM.write_report(template, ... % template file
                    report_name, default.show_report, 0); % output format
            else
                ksCOM.report(default.show_report);
            end
        end
        
        function SH = Pitting(obj, varargin)
            % Safety factor for surface durability (against pitting):
            SH = obj.S_Hmin*obj.sigma_HP./obj.sigma_H; % pinion/planet
            
            if(isrow(SH))
                SH = SH';
            end
        end
        
    end
    
    methods(Static)
        function calc = example_01()
            %EXAMPLE_01 runs the first calculation example from 
            %                      Parameter                         Symbol         Problem       Reference         Calculated         Rel_Error_pc
            %     ___________________________________________    _______________    _______    _______________    _______________    ________________
            %     'Minimum safety factor (pitting)'              'S_Hmin'           'no'       [  1.0000e+000]    [  1.0000e+000]    [   0.0000e+000]
            %     'Pitting safety factor for pinion'             'S_H1'             'YES'      [  1.0285e+000]    [  1.1102e+000]    [  -7.9383e+000]
            %     'Pitting safety factor for wheel'              'S_H2'             'YES'      [  1.0870e+000]    [  1.1733e+000]    [  -7.9388e+000]
            %     'Permissible contact stress for pinion'        'sigma_HP1'        'YES'      [  1.3385e+003]    [  1.3626e+003]    [  -1.8010e+000]
            %     'Permissible contact stress for wheel'         'sigma_HP2'        'YES'      [  1.4145e+003]    [  1.4400e+003]    [  -1.8010e+000]
            %     'Nominal contact stress'                       'sigma_H0'         'no'       [  1.2066e+003]    [  1.2070e+003]    [ -37.1763e-003]
            %     'Contact stress for pinion'                    'sigma_H1'         'YES'      [  1.3014e+003]    [  0.0000e+000]    [ 100.0000e+000]
            %     'Contact stress for wheel'                     'sigma_H2'         'YES'      [  1.3014e+003]    [  0.0000e+000]    [ 100.0000e+000]
            %     'Nominal tangential load, [N]'                 'F_t'              'no'       [127.3520e+003]    [127.3524e+003]    [-299.5607e-006]
            %     'Pitch line velocity, [m/s]'                   'v'                'no'       [  2.6640e+000]    [  2.6642e+000]    [  -7.4461e-003]
            %     'Dynamic factor'                               'K_v'              'no'       [  1.0030e+000]    [  1.0055e+000]    [-249.2935e-003]
            %     'Transverse load factor (root stress)'         'K_Falpha'         'no'       [  1.0000e+000]    [  1.0000e+000]    [   0.0000e+000]
            %     'Transverse load factor (contact stress)'      'K_Halpha'         'no'       [  1.0000e+000]    [  1.0000e+000]    [   0.0000e+000]
            %     'Face load factor (root stress)'               'K_Fbeta'          'YES'      [  1.1280e+000]    [  1.0229e+000]    [   9.3189e+000]
            %     'Face load factor (contact stress)'            'K_Hbeta'          'YES'      [  1.1600e+000]    [  1.0283e+000]    [  11.3525e+000]
            %     'Virtual number of teeth for pinion'           'z_n1'             'no'       [ 18.9050e+000]    [ 18.9051e+000]    [-650.6211e-006]
            %     'Virtual number of teeth for wheel'            'z_n2'             'no'       [114.5430e+000]    [114.5428e+000]    [ 171.0638e-006]
            %     'Number of load cycles for pinion'             'N_L1'             'no'       [  1.0800e+009]    [  1.0800e+009]    [   0.0000e+000]
            %     'Number of load cycles for wheel'              'N_L2'             'no'       [178.3000e+006]    [178.2524e+006]    [  26.6813e-003]
            %     'Life factor for pinion'                       'Z_NT1'            'no'       [910.0000e-003]    [910.0545e-003]    [  -5.9848e-003]
            %     'Life factor for wheel'                        'Z_NT2'            'no'       [962.0000e-003]    [961.7587e-003]    [  25.0846e-003]
            %     'Size factor'                                  'Z_X'              'no'       [  1.0000e+000]    [  1.0000e+000]    [   0.0000e+000]
            %     'Work hardening factor'                        'Z_W'              'no'       [  1.0000e+000]    [  1.0000e+000]    [   0.0000e+000]
            %     'Roughness factor'                             'Z_R'              'YES'      [965.9900e-003]    [983.3849e-003]    [  -1.8007e+000]
            %     'Velocity factor'                              'Z_v'              'no'       [969.1100e-003]    [969.1079e-003]    [ 216.3428e-006]
            %     'Lubricant factor'                             'Z_L'              'no'       [  1.0474e+000]    [  1.0474e+000]    [-350.1084e-006]
            %     'Helix angle factor'                           'Z_beta'           'no'       [  1.0194e+000]    [  1.0194e+000]    [-366.7700e-006]
            %     'Contact ratio factor'                         'Z_eps'            'no'       [803.0000e-003]    [803.7539e-003]    [ -93.8840e-003]
            %     'Elasticity factor'                            'Z_E'              'no'       [189.8117e+000]    [189.8117e+000]    [-210.7352e-009]
            %     'Single pair tooth contact factor'             'Z_B'              'no'       [  1.0000e+000]    [  1.0000e+000]    [   0.0000e+000]
            %     'Single pair tooth contact factor'             'Z_D'              'no'       [  1.0000e+000]    [  1.0000e+000]    [   0.0000e+000]
            %     'Zone factor'                                  'Z_H'              'no'       [  2.3953e+000]    [  2.3953e+000]    [-179.4742e-006]
            %     'Single stiffness, [N/(mm-um)]'                'c''               'no'       [ 12.3705e+000]    [ 12.3620e+000]    [  68.2642e-003]
            %     'Theoretical single stiffness, [N/(mm-um)]'    'c_th''            'no'       [ 17.8558e+000]    [ 17.8436e+000]    [  68.2770e-003]
            %     'Mesh stiffness, [N/(mm-um)]'                  'c_gamma_alpha'    'no'       [ 17.4648e+000]    [ 17.4553e+000]    [  54.8918e-003]
            %     'Mesh stiffness, [N/(mm-um)]'                  'c_gamma_beta'     'no'       [ 14.8451e+000]    [ 14.8370e+000]    [  54.8750e-003]
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
            
            calc = ISO_6336_KS(gset, 'calculation', 'default', ...
                                     'P_rated'    , P_r, ...
                                     'S_Hmin'     , 1.0, ...
                                     'S_Fmin'     , 1.0, ...
                                     'L_h'        , 50.0e3, ...
                                     'K_A'        , 1.0, ...
                                     'n_nominal'  , [n_1 NaN], ...
                                     'lubricant_ID', 11180, ...
                                     'C_a'         , 70.0);

            tab_set = {'Minimum safety factor (pitting)'          , 'S_Hmin'       , calc.S_Hmin       , 1.0;
                       'Pitting safety factor for pinion'         , 'S_H1'         , calc.S_H(1)       , 1.02853;
                       'Pitting safety factor for wheel'          , 'S_H2'         , calc.S_H(2)       , 1.08696;
                       'Permissible contact stress for pinion'    , 'sigma_HP1'    , calc.sigma_HP(1)  , 1338.48050;
                       'Permissible contact stress for wheel'     , 'sigma_HP2'    , calc.sigma_HP(2)  , 1414.52551;
                       'Nominal contact stress'                   , 'sigma_H0'     , calc.sigma_H0(1)  , 1206.58207;
                       'Contact stress for pinion'                , 'sigma_H1'     , calc.sigma_H(1)   , 1301.35343;
                       'Contact stress for wheel'                 , 'sigma_H2'     , calc.sigma_H(2)   , 1301.35343;
                       'Nominal tangential load, [N]'             , 'F_t'          , calc.F_t          , 127352.0;
                       'Pitch line velocity, [m/s]'               , 'v'            , calc.v_pitch_line , 2.664;
                       'Dynamic factor'                           , 'K_v'          , calc.K_v          , 1.003;
                       'Transverse load factor (root stress)'     , 'K_Falpha'     , calc.K_Falpha(1)  , 1.0;
                       'Transverse load factor (contact stress)'  , 'K_Halpha'     , calc.K_Halpha(1)  , 1.0;
                       'Face load factor (root stress)'           , 'K_Fbeta'      , calc.K_Fbeta(1)   , 1.12803;
                       'Face load factor (contact stress)'        , 'K_Hbeta'      , calc.K_Hbeta(1)   , 1.16;
                       'Virtual number of teeth for pinion'       , 'z_n1'         , calc.z_n(1)       , 18.905;
                       'Virtual number of teeth for wheel'        , 'z_n2'         , calc.z_n(2)       , 114.543;
                       'Number of load cycles for pinion'         , 'N_L1'         , calc.N_L(1)       , 1.080e9;
                       'Number of load cycles for wheel'          , 'N_L2'         , calc.N_L(2)       , 1.783e8;
                       'Life factor for pinion'                   , 'Z_NT1'        , calc.Z_NT(1)      , 0.91;
                       'Life factor for wheel'                    , 'Z_NT2'        , calc.Z_NT(2)      , 0.962;
                       'Size factor'                              , 'Z_X'          , calc.Z_X          , 1.0;
                       'Work hardening factor'                    , 'Z_W'          , calc.Z_W          , 1.0;
                       'Roughness factor'                         , 'Z_R'          , calc.Z_R(1)       , 0.96599;
                       'Velocity factor'                          , 'Z_v'          , calc.Z_v          , 0.96911;
                       'Lubricant factor'                         , 'Z_L'          , calc.Z_L          , 1.04739;
                       'Helix angle factor'                       , 'Z_beta'       , calc.Z_beta(1)    , 1.01944;
                       'Contact ratio factor'                     , 'Z_eps'        , calc.Z_eps(1)     , 0.803;
                       'Elasticity factor'                        , 'Z_E'          , calc.Z_E(1)       , 189.81170;
                       'Single pair tooth contact factor'         , 'Z_B'          , calc.Z_BD(1)      , 1.0;
                       'Single pair tooth contact factor'         , 'Z_D'          , calc.Z_BD(2)      , 1.0;
                       'Zone factor'                              , 'Z_H'          , calc.Z_H(1)       , 2.39533;
                       'Single stiffness, [N/(mm-um)]'            , 'c'''          , calc.cprime       , 12.37047;
                       'Theoretical single stiffness, [N/(mm-um)]', 'c_th'''       , calc.cprime_th    , 17.85584;
                       'Mesh stiffness, [N/(mm-um)]'              , 'c_gamma_alpha', calc.c_gamma_alpha, 17.46485;
                       'Mesh stiffness, [N/(mm-um)]'              , 'c_gamma_beta' , calc.c_gamma_beta , 14.84512;
                       };

            Parameter = tab_set(:, 1);
            Symbol    = tab_set(:, 2);
            Calculated  = tab_set(:, 3);
            Reference   = tab_set(:, 4);
            Rel_Error_pc   = (1.0 - cell2mat(Calculated)./cell2mat(Reference))*100.0;
            Problem = repmat("no", length(Rel_Error_pc), 1);
            Problem(abs(Rel_Error_pc) > 1.0) = "YES";
            Problem = cellstr(Problem);
            Rel_Error_pc = mat2cell(Rel_Error_pc, ones(length(Rel_Error_pc),1));
            
            tab = table(Parameter, Symbol, Problem, Reference, Calculated, Rel_Error_pc);
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
            
            calc = ISO_6336(gset, 'calculation', 'default', ...
                                  'P_rated'    , P_r      , ...
                                  'S_Hmin'     , 1.0      , ...
                                  'S_Fmin'     , 1.0         , ...
                                  'L_h'        , 1.0e3, ...
                                  'K_A'        , 1.0      , ...
                                  'n_nominal'  , [n_1 NaN], ...
                                  'nu_40'      , 320.0    , ...
                                  'C_a'        , 70.0);

        end
        
    end

    %% Set methods:
    methods
    end
    
    %% Get methods:
    methods
        function val = get.S_H(obj)
            n = length(obj.z);
            val = zeros(1, n);
            
            for idx = 1:n
                val(idx) = obj.ks.get_var(sprintf('ZPP[%d].Flanke.SH', idx - 1));
            end
        end
        
        function val = get.S_F(obj)
            n = length(obj.z);
            val = zeros(1, n);
            
            for idx = 1:n
                val(idx) = obj.ks.get_var(sprintf('ZPP[%d].Fuss.SF', idx - 1));
            end
        end
        
        function val = get.sigma_H0(obj)
            % Nominal contact stress at pitch point:
            n = length(obj.z);
            val = zeros(1, n);
            
            for idx = 1:n
                val(idx) = obj.ks.get_var(sprintf('ZP[%d].Flanke.sigH0', idx - 1));
            end
        end
        
        function val = get.sigma_H(obj)
            n = length(obj.z);
            val = zeros(1, n);

            for idx = 1:n
%                 val(idx) = obj.ks.get_var(sprintf('ZP[%d].flankbreak.resultsForDocumentation.sigmaH', idx - 1));
                val(idx) = obj.ks.get_var(sprintf('ZP[%d].Flanke.sigH', idx - 1));
            end
        end
        
        function val = get.sigma_HP(obj)
            % Permissible contact stress:
            n = length(obj.z);
            val = zeros(1, n);

            for idx = 1:n
                % that's actually the pitting stress limit, sigma_HG:
                val(idx) = obj.ks.get_var(sprintf('ZPP[%d].Flanke.sigHP', idx - 1));
            end
            val = val./obj.ks.get_var('ZS.SSi.Flanke');
        end
        
        function val = get.sigma_Hlim(obj)
            n = length(obj.z);
            val = zeros(1, n);

            for idx = 1:n
                val(idx) = obj.ks.get_var(sprintf('ZR[%d].mat.limh', idx - 1));
            end
        end
        
        function val = get.F_t(obj)
            % [N], Nominal tangential load: ZPP[0].Schn.Ft
            n = length(obj.z);
            val = zeros(1, n);

            for idx = 1:n
%                 val(idx) = obj.ks.get_var(sprintf('ZPP[%d].Schn.Ft', idx - 1));
                val(idx) = obj.ks.get_var(sprintf('ZP[%d].Ft', idx - 1));
            end
        end
        
        function val = get.Z_NT(obj)
            n = length(obj.z);
            val = zeros(1, n);
            
            for idx = 1:n
                val(idx) = obj.ks.get_var(sprintf('ZR[%d].Flanke.ZNT', idx - 1));
            end
        end
        
        function val = get.N_L(obj)
            n = length(obj.z);
            val = zeros(1, n);
            
            for idx = 1:n
                val(idx) = obj.ks.get_var(sprintf('ZR[%d].NL', idx - 1));
            end
        end
        
        function val = get.K_v(obj)
            val = obj.ks.get_var('ZS.KVcalc');
            % or the max of obj.ks.get_var(sprintf('ZP[%d].KV.KV', idx - 1))
        end
        
        function val = get.line_load(obj)
            val = obj.F_t.*obj.K_A.*obj.K_gamma./obj.b;
        end
        
        function val = get.v_pitch_line(obj)
            val = obj.ks.get_var('ZS.v');
        end
        
        function val = get.K_Halpha(obj)
            n = length(obj.z);
            val = zeros(1, n);

            for idx = 1:n
                val(idx) = obj.ks.get_var(sprintf('ZP[%d].KHa', idx - 1));
            end
        end
        
        function val = get.K_Falpha(obj)
            n = length(obj.z);
            val = zeros(1, n);

            for idx = 1:n
                val(idx) = obj.ks.get_var(sprintf('ZP[%d].KFa', idx - 1));
            end
        end
        
        function val = get.K_Hbeta(obj)
            n = length(obj.z);
            val = zeros(1, n);

            for idx = 1:n
                val(idx) = obj.ks.get_var(sprintf('ZP[%d].KHb', idx - 1));
            end
        end
        
        function val = get.K_Fbeta(obj)
            n = length(obj.z);
            val = zeros(1, n);

            for idx = 1:n
                val(idx) = obj.ks.get_var(sprintf('ZP[%d].KFb', idx - 1));
            end
        end
        
        function val = get.Z_eps(obj)
            n = length(obj.z);
            val = zeros(1, n);

            for idx = 1:n
                val(idx) = obj.ks.get_var(sprintf('ZP[%d].Flanke.Zeps', idx - 1));
            end
        end
        
        function val = get.Z_beta(obj)
            n = length(obj.z);
            val = zeros(1, n);

            for idx = 1:n
                val(idx) = obj.ks.get_var(sprintf('ZP[%d].Flanke.Zbet', idx - 1));
            end
        end
        
        function val = get.Z_BD(obj)
            n = length(obj.z);
            val = zeros(1, n);

            for idx = 1:n
                val(idx) = obj.ks.get_var(sprintf('ZPP[%d].Flanke.ZBD', idx - 1));
            end
        end
        
        function val = get.Z_E(obj)
            n = length(obj.z);
            val = zeros(1, n);

            for idx = 1:n
                val(idx) = obj.ks.get_var(sprintf('ZP[%d].Flanke.ZE', idx - 1));
            end
        end
        
        function val = get.Z_H(obj)
            n = length(obj.z);
            val = zeros(1, n);

            for idx = 1:n
                val(idx) = obj.ks.get_var(sprintf('ZP[%d].Flanke.ZH', idx - 1));
            end
        end
        
        function val = get.Z_R(obj)
            n = length(obj.z);
            val = zeros(1, n);

            for idx = 1:n
                val(idx) = obj.ks.get_var(sprintf('ZPP[%d].Flanke.ZR', idx - 1));
            end
        end
        
        function val = get.Z_L(obj)
            n = length(obj.z);
            val = zeros(1, n);

            for idx = 1:n
                val(idx) = obj.ks.get_var(sprintf('ZPP[%d].Flanke.ZL', idx - 1));
            end
        end
        
        function val = get.Z_v(obj)
            n = length(obj.z);
            val = zeros(1, n);

            for idx = 1:n
                val(idx) = obj.ks.get_var(sprintf('ZPP[%d].Flanke.ZV', idx - 1));
            end
        end
        
        function val = get.Z_W(obj)
            n = length(obj.z);
            val = zeros(1, n);

            for idx = 1:n
                val(idx) = obj.ks.get_var(sprintf('ZPP[%d].Flanke.ZW', idx - 1));
            end
        end
        
        function val = get.Z_X(obj)
            n = length(obj.z);
            val = zeros(1, n);

            for idx = 1:n
                val(idx) = obj.ks.get_var(sprintf('ZPP[%d].Flanke.ZX', idx - 1));
            end
        end
        
    end
end
