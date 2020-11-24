classdef MATISO_6336 < ISO_6336
    %MATISO_6336 do some calculations following the ISO 6336 standard.
    % see also: ISO_6336, KSISO_6336, Gear_Set
    %
    % References:
    %   [1] Calculation of load capacity of spur and helical gears -- 
    % Part 1: Basic principles, introduction and general influence factors,
    % ISO 6336-1:2019, 2019.
    %   [2] Calculation of load capacity of spur and helical gears -- 
    % Part 2: Calculation of surface durability (pitting), 
    % ISO 6336-2:2019, 2019. 
    %   [3] Calculation of load capacity of spur and helical gears -- 
    % Part 6: Calculation of service life under variable load, 
    % ISO 6336-6:2019, 2019.
    %   [4] Calculation of load capacity of spur and helical gears -- 
    % Calculation examples for the application of ISO 6336 parts 1, 2, 3, 5,
    % ISO/TR 6336-30:2017, 2017.
    %   [5] Arnaudov, K., Karaivanov, D. (2019). Planetary Gear Trains. 
    % Boca Raton: CRC Press, https://doi.org/10.1201/9780429458521
    %
    
    properties(Access = public)
        nu_40;       % [mm^2/s], Kinematic viscosity of the lubricant at 40 deg. C
        C_a;         % [um],     Tip relief
    end
    
    properties(Access = private)
        K_v_val;
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
        function obj = MATISO_6336(gset, varargin)
            if(~exist('gset', "var"))
                gset = Gear_Set();
            elseif(~strcmp(class(gset), "Gear_Set"))
                error('gset is not a Gear_Set object.');
            end
            
            T_1 = 9.0e3;
            n_1 = 360.0;
            P_r = 1.0e-3*T_1*n_1*pi/30.0;
            
            C_ay = (1.0/18.0)*([gset.material(1:2).sigma_Hlim]./97.0 - 18.45).^2 + 1.5;
            C_ay = mean(C_ay);

            default = {'P_rated'    , P_r, ...
                       'S_Hmin'     , 1.0, ...
                       'S_Fmin'     , 1.0, ...
                       'L_h'        , 50.0e3, ...
                       'K_A'        , 1.0, ...
                       'nu_40'      , 320.0, ...
                       'n_nominal'  , [n_1, NaN], ...
                       'C_a'        , C_ay};
            
            default = scaling_factor.process_varargin(default, varargin);
            
            obj@ISO_6336(gset, 'P_rated'    , default.P_rated, ...
                               'S_Hmin'     , default.S_Hmin, ...
                               'S_Fmin'     , default.S_Fmin, ...
                               'L_h'        , default.L_h, ...
                               'K_A'        , default.K_A, ...
                               'n_nominal'  , default.n_nominal);
            
            obj.nu_40   = default.nu_40;
            obj.C_a     = default.C_a;

            obj.K_v_val = obj.dynamic_factor();
        end
        
        function save_KS(obj, name)
            % do nothing...
        end
    end
    
    %% Calculation methods:
    methods
        function [SH, SF] = safety_factors(obj, varargin)
            SH = obj.Pitting(varargin{:});
            SF = nan(size(SH));
        end
        
        function SH = Pitting(obj, varargin)
            % Safety factor for surface durability (against pitting):
            SH = obj.S_Hmin*obj.sigma_HP./obj.sigma_H; % pinion/planet
            
            if(isrow(SH))
                SH = SH';
            end
        end
        
%         function SF = tooth_bending(obj, varargin)
%             % to do...
%         end

        function ZN = life_factor(obj, N)
            
            ZN = zeros(2, length(N));
            
            for idx = 1:2
                line = obj.material(idx).row;
                ratio = obj.Z_NT(idx); % obj.sigma_HP_stat(idx)/obj.sigma_HP_ref(idx);
                
                switch line
                    case 1
                        % St, V, GGG (perl. bai.), GTS (perl.), Eh, IF (when limited pitting is permitted)
                        val = power(3.0e8./6.0e5, 0.3705*log(ratio)).*(               (N <  6.0e5)) + ... % if 1
                              power(3.0e8./N    , 0.3705*log(ratio)).*((6.0e5 <  N) & (N <= 1.0e7)) + ... % if 2
                              power(1.0e9./N    , 0.2791*log(ratio)).*((1.0e7 <  N) & (N <= 1.0e9)) + ...
                              1.0                                   .*((1.0e9 <  N));
                        
                    case 2
                        % St, V, GGG (perl. bai.), GTS (perl.), Eh, IF
                        val = power(5.0e7./1.0e5, 0.3705*log(ratio)).*(               (N <  1.0e5)) + ... % if 1
                              power(5.0e7./N    , 0.3705*log(ratio)).*((1.0e5 <  N) & (N <= 5.0e7)) + ... % if 2
                              1.0                                   .*((5.0e7 <  N));
                        
                    case {3, 4}
                        % GG, GGG (ferr.), NT (nitr.), NV (nitr.)
                        % NV (nitrocar.)
                        val = power(2.0e6./1.0e5, 0.7686*log(ratio)).*(               (N <  1.0e5)) + ... % if 1
                              power(2.0e6./N    , 0.7686*log(ratio)).*((1.0e5 <  N) & (N <= 2.0e6)) + ... % if 2
                              1.0                                   .*((2.0e6 <  N));
                        
                    otherwise
                        error("MATISO_6336:ZN", "Invalid input [%d].\nValid options are 1 to 4.\n", line);
                end
                ZN(idx, :) = val;
            end
            
        end
        
    end
    
    methods(Static)
        function calc = example_01()
            %EXAMPLE_01 runs the worked example 01 from [4], entitled 
            % 'Single helical case carburized gear pair', and compares the 
            % results against the values found in [4], Annex A.
            %
            % Results from 13.11.2020:
            % - Z_NT was obtained using y_lim = 0.85, which is not the case
            % for other systems.
            %                       Parameter                           Symbol          Problem        Reference           Calculated           Rel_Error_pc   
            %     _____________________________________________    _________________    _______    _________________    _________________    __________________
            % 
            %     {'Minimum safety factor (pitting)'          }    {'S_Hmin'       }      ""       {[  1.0000e+000]}    {[  1.0000e+000]}    {[   0.0000e+000]}
            %     {'Pitting safety factor for pinion'         }    {'S_H1'         }      ""       {[  1.0285e+000]}    {[  1.0299e+000]}    {[-137.6429e-003]}
            %     {'Pitting safety factor for wheel'          }    {'S_H2'         }      ""       {[  1.0870e+000]}    {[  1.0885e+000]}    {[-138.1352e-003]}
            %     {'Permissible contact stress for pinion'    }    {'sigma_HP1'    }      ""       {[  1.3385e+003]}    {[  1.3385e+003]}    {[ -30.2347e-009]}
            %     {'Permissible contact stress for wheel'     }    {'sigma_HP2'    }      ""       {[  1.4145e+003]}    {[  1.4145e+003]}    {[ 328.5801e-009]}
            %     {'Nominal contact stress'                   }    {'sigma_H0'     }      ""       {[  1.2066e+003]}    {[  1.2065e+003]}    {[   8.1404e-003]}
            %     {'Contact stress for pinion'                }    {'sigma_H1'     }      ""       {[  1.3014e+003]}    {[  1.2996e+003]}    {[ 137.4942e-003]}
            %     {'Contact stress for wheel'                 }    {'sigma_H2'     }      ""       {[  1.3014e+003]}    {[  1.2996e+003]}    {[ 137.4942e-003]}
            %     {'Nominal tangential load, [N]'             }    {'F_t'          }      ""       {[127.3520e+003]}    {[127.3524e+003]}    {[-299.5607e-006]}
            %     {'Pitch line velocity, [m/s]'               }    {'v'            }      ""       {[  2.6640e+000]}    {[  2.6642e+000]}    {[  -7.4461e-003]}
            %     {'Dynamic factor'                           }    {'K_v'          }      ""       {[  1.0030e+000]}    {[  1.0005e+000]}    {[ 246.6092e-003]}
            %     {'Transverse load factor (root stress)'     }    {'K_Falpha'     }      ""       {[  1.0000e+000]}    {[  1.0000e+000]}    {[   0.0000e+000]}
            %     {'Transverse load factor (contact stress)'  }    {'K_Halpha'     }      ""       {[  1.0000e+000]}    {[  1.0000e+000]}    {[   0.0000e+000]}
            %     {'Face load factor (root stress)'           }    {'K_Fbeta'      }      ""       {[  1.1280e+000]}    {[  1.1281e+000]}    {[  -5.4065e-003]}
            %     {'Face load factor (contact stress)'        }    {'K_Hbeta'      }      ""       {[  1.1600e+000]}    {[  1.1596e+000]}    {[  30.8898e-003]}
            %     {'Virtual number of teeth for pinion'       }    {'z_n1'         }      ""       {[ 18.9050e+000]}    {[ 18.9051e+000]}    {[-650.6211e-006]}
            %     {'Virtual number of teeth for wheel'        }    {'z_n2'         }      ""       {[114.5430e+000]}    {[114.5428e+000]}    {[ 171.0638e-006]}
            %     {'Number of load cycles for pinion'         }    {'N_L1'         }      ""       {[  1.0800e+009]}    {[  1.0800e+009]}    {[   0.0000e+000]}
            %     {'Number of load cycles for wheel'          }    {'N_L2'         }      ""       {[178.3000e+006]}    {[178.2524e+006]}    {[  26.6813e-003]}
            %     {'Life factor for pinion'                   }    {'Z_NT1'        }      ""       {[910.0000e-003]}    {[910.0545e-003]}    {[  -5.9848e-003]}
            %     {'Life factor for wheel'                    }    {'Z_NT2'        }      ""       {[962.0000e-003]}    {[961.7587e-003]}    {[  25.0846e-003]}
            %     {'Size factor'                              }    {'Z_X'          }      ""       {[  1.0000e+000]}    {[  1.0000e+000]}    {[   0.0000e+000]}
            %     {'Work hardening factor'                    }    {'Z_W'          }      ""       {[  1.0000e+000]}    {[  1.0000e+000]}    {[   0.0000e+000]}
            %     {'Roughness factor'                         }    {'Z_R'          }      ""       {[965.9900e-003]}    {[965.9878e-003]}    {[ 232.3318e-006]}
            %     {'Velocity factor'                          }    {'Z_v'          }      ""       {[969.1100e-003]}    {[969.1142e-003]}    {[-434.0740e-006]}
            %     {'Lubricant factor'                         }    {'Z_L'          }      ""       {[  1.0474e+000]}    {[  1.0474e+000]}    {[ 368.0474e-006]}
            %     {'Helix angle factor'                       }    {'Z_beta'       }      ""       {[  1.0194e+000]}    {[  1.0194e+000]}    {[-366.7828e-006]}
            %     {'Contact ratio factor'                     }    {'Z_eps'        }      ""       {[803.0000e-003]}    {[803.3898e-003]}    {[ -48.5416e-003]}
            %     {'Elasticity factor'                        }    {'Z_E'          }      ""       {[189.8117e+000]}    {[189.8117e+000]}    {[-230.5266e-009]}
            %     {'Single pair tooth contact factor'         }    {'Z_B'          }      ""       {[  1.0000e+000]}    {[  1.0000e+000]}    {[   0.0000e+000]}
            %     {'Single pair tooth contact factor'         }    {'Z_D'          }      ""       {[  1.0000e+000]}    {[  1.0000e+000]}    {[   0.0000e+000]}
            %     {'Zone factor'                              }    {'Z_H'          }      ""       {[  2.3953e+000]}    {[  2.3953e+000]}    {[-179.4571e-006]}
            %     {'Single stiffness, [N/(mm-um)]'            }    {'c''           }      ""       {[ 12.3705e+000]}    {[ 12.3620e+000]}    {[  68.2642e-003]}
            %     {'Theoretical single stiffness, [N/(mm-um)]'}    {'c_th''        }      ""       {[ 17.8558e+000]}    {[ 17.8436e+000]}    {[  68.2770e-003]}
            %     {'Mesh stiffness, [N/(mm-um)]'              }    {'c_gamma_alpha'}      ""       {[ 17.4648e+000]}    {[ 17.4553e+000]}    {[  54.8918e-003]}
            %     {'Mesh stiffness, [N/(mm-um)]'              }    {'c_gamma_beta' }      ""       {[ 14.8451e+000]}    {[ 14.8370e+000]}    {[  54.8750e-003]}
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
            
            calc = MATISO_6336(gset, 'P_rated'    , P_r, ...
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
        
        function calc = example_load_spectra()
            %EXAMPLE_LOAD_SPECTRA runs the example in [3], Annex C,
            % entitled 'Example calculation for safety factor from given
            % load spectrum', and compares the results against the values
            % found in [3], Annex A.
            %
            
%
%                        'Dynamic factor'                           , 'K_v'          , calc.K_v          , 1.003;
%                        'Nominal tangential load, [N]'             , 'F_t'          , calc.F_t          , 127352.0;
%                        'Pitch line velocity, [m/s]'               , 'v'            , calc.v_pitch_line , 2.664;
%                        'Transverse load factor (root stress)'     , 'K_Falpha'     , calc.K_Falpha     , 1.0;
%                        'Transverse load factor (contact stress)'  , 'K_Halpha'     , calc.K_Halpha     , 1.0;
%                        'Face load factor (root stress)'           , 'K_Fbeta'      , calc.K_Fbeta      , 1.12803;
%                        'Face load factor (contact stress)'        , 'K_Hbeta'      , calc.K_Hbeta      , 1.16;
%                        'Virtual number of teeth for pinion'       , 'z_n1'         , calc.z_n(1)       , 18.905;
%                        'Virtual number of teeth for wheel'        , 'z_n2'         , calc.z_n(2)       , 114.543;
%                        'Number of load cycles for pinion'         , 'N_L1'         , calc.N_L(1)       , 1.080e9;
%                        'Number of load cycles for wheel'          , 'N_L2'         , calc.N_L(2)       , 1.783e8;
%                        'Helix angle factor'                       , 'Z_beta'       , calc.Z_beta       , 1.01944;
%                        'Contact ratio factor'                     , 'Z_eps'        , calc.Z_eps        , 0.803;
%                        'Elasticity factor'                        , 'Z_E'          , calc.Z_E          , 189.81170;
%                        'Single pair tooth contact factor'         , 'Z_B'          , calc.Z_BD(1)      , 1.0;
%                        'Single pair tooth contact factor'         , 'Z_D'          , calc.Z_BD(2)      , 1.0;
%                        'Zone factor'                              , 'Z_H'          , calc.Z_H          , 2.39533;
%                        'Single stiffness, [N/(mm-um)]'            , 'c'''          , calc.cprime       , 12.37047;
%                        'Theoretical single stiffness, [N/(mm-um)]', 'c_th'''       , calc.cprime_th    , 17.85584;
%                        'Mesh stiffness, [N/(mm-um)]'              , 'c_gamma_alpha', calc.c_gamma_alpha, 17.46485;
%                        'Mesh stiffness, [N/(mm-um)]'              , 'c_gamma_beta' , calc.c_gamma_beta , 14.84512;
%

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

            T_1 = 25423.0; % [N-m], from [6], Table 4, bin 3
            n_1 = 35.2;    % [1.0/min.]
            P_r = 1.0e-3*T_1*n_1*pi/30.0;
            
            calc = MATISO_6336(gset, 'P_rated'    , P_r      , ...
                                     'S_Hmin'     , 1.0      , ...
                                     'S_Fmin'     , 1.0      , ...
                                     'L_h'        , 1.0e3    , ...
                                     'K_A'        , 1.0      , ...
                                     'n_nominal'  , [n_1 NaN], ...
                                     'C_a'        , 70.0);
            
            tab_set = {'Minimum safety factor (pitting)'          , 'S_Hmin'       , calc.S_Hmin       , 1.0;
                       'Pitting safety factor for pinion'         , 'S_H1'         , calc.S_H(1)       , 1.374;
                       'Permissible contact stress for pinion'    , 'sigma_HP1'    , calc.sigma_HP(1)  , 2400.0;
                       'Nominal contact stress'                   , 'sigma_H0'     , calc.sigma_H0     , 1206.58207;
                       'Contact stress for pinion'                , 'sigma_H1'     , calc.sigma_H(1)   , 1747.0;
                       'Life factor for pinion'                   , 'Z_NT1'        , calc.Z_NT(1)      , 0.91;
                       'Size factor'                              , 'Z_X'          , calc.Z_X          , 1.0;
                       'Work hardening factor'                    , 'Z_W'          , calc.Z_W          , 1.0;
                       'Roughness factor'                         , 'Z_R'          , calc.Z_R          , 1.01;
                       'Velocity factor'                          , 'Z_v'          , calc.Z_v          , 0.943;
                       'Lubricant factor'                         , 'Z_L'          , calc.Z_L          , 1.02;
                       };

            tab_04 = readtable('table_04.csv');
            torque_edges = tab_04.Torque_Max_Nm;
            speed_Max = n_1*ones(size(torque_edges));
            Max_power = torque_edges.*1.0e-3 .* ... % [kN-m]
                        speed_Max.*(pi/30.0);     % [rad/s] = [kW]

            L_h30y = 30.0*365.25*24.0;

            n_30y = tab_04.load_cycles*L_h30y/sum(tab_04.time_h);

            calc_tmp = calc;
            calc_tmp.K_A = 1.0;

            Khb    = zeros(size(speed_Max));
            SH     = zeros(size(speed_Max));
            ZNT    = zeros(size(speed_Max));
            sigmaH = zeros(2, size(speed_Max, 1));

            range = 3:44; % Non-zero values
            for idx = range %1:length(speed_Max)
                calc_tmp.P_rated = Max_power(idx);
                calc_tmp.n_nominal = calc_tmp.speed_ratio*abs(speed_Max(idx));
                calc_tmp.L_h = sum(n_30y(1:idx));

                sigmaH(:, idx) = calc_tmp.sigma_H';
                Khb(idx)       = calc_tmp.K_Hbeta;
                SH(idx)        = calc_tmp.S_H(1);
                ZNT(idx)       = calc_tmp.Z_NT(1);
            end
            
            sigmaHG = sigmaH.*1.275;
            
%             Ni = nan(size(SH));
%             Ui = Ni;

%             tab_C2_calc = table(tab_04.num_bin, torque_edges*1.0e-3, tab_04.time_s       , speed_Max         , n_30y       , Khb      , sigmaH.*SH   , ZNT   , Ni   , Ui, ...
%               'variableNames', {'Bin_No'      , 'Pinion_Torque'    , 'Operating_time_70d', 'Pinion_speed_rpm', 'Cycles_30y', 'K_Hbeta', 'sigma_HxS_H', 'Z_NT', 'N_i', 'U_i'});

            num_cycle = [0.0; cumsum(n_30y)]';

            [NN, sigHH] = calc.Pitting_SN_curve();

            %% Plotting LDD:
            torque_hist = [torque_edges; tab_04.Torque_min_Nm(end)];
            torque_hist = flipud(torque_hist)*1.0e-3;
            zdx = num_cycle == 0.0;
            cycle_num = num_cycle;
            cycle_num(zdx) = 1;

            figure('unit', 'centimeters', 'position', [5, 5, 23, 16]);
            subplot(211)
            histogram('binEdges' , torque_hist', ...
                      'binCounts', flipud(tab_04.load_cycles));
            title('ISO 6336-6:2019, Table 4 -- Histogram')
            xlabel('T, [kN-m]');
            ylabel('N, [-]');

            subplot(212)
            histogram('binEdges' , cycle_num, ...
                      'binCounts', torque_edges*1.0e-3);
            title('ISO 6336-6:2019, Table 4 -- LDD');
            set(gca, 'xScale', 'log')
            xlabel('N, [-]');
            ylabel('T, [kN-m]');
            xlim([1e2 1e6]);

            fig_axes = findobj(gcf, 'Type', 'Axes');
            set(fig_axes, 'fontName', 'Times');
            set(fig_axes, 'fontSize', 12);
            set(fig_axes, 'box', 'on');

            %% Figure C.1: Load Spectra and S-N Curve
            fig_C1 = readtable('data_Fig_C1.csv');
            
            figure('unit', 'centimeters', 'position', [5, 5, 23, 16]);
            subplot(211)
            hold on;
            plot(NN, sigHH(1, :), 'k', 'lineWidth', 2.0);
            histogram('binEdges' , cycle_num, ...
                      'binCounts', sigmaH(1, :)*1.275, ...
                      'displayStyle', 'stairs', ...
                      'edgeColor', 'r', ...
                      'lineWidth', 2.0);
            histogram('binEdges' , cycle_num, ...
                      'binCounts', sigmaH(2, :), ...
                      'displayStyle', 'stairs', ...
                      'edgeColor', 'b', ...
                      'lineWidth', 2.0);
            title('ISO 6336-6:2019, Fig. C.1 -- Load Spectra and S-N Curve (Obtained)')
            ylabel('\sigma_H, [N/mm^2]');

            subplot(212)
            hold on;
            plot(fig_C1.N_SN, fig_C1.sigma_SN, 'k-', 'lineWidth', 2.0);
            plot(fig_C1.N_H1, fig_C1.sigma_H1, 'r', 'lineWidth', 2.0);
            plot(fig_C1.N_H2, fig_C1.sigma_H2, 'b', 'lineWidth', 2.0);
            xlabel('N, [-]');
            ylabel('\sigma_H, [N/mm^2]');
            title('ISO 6336-6:2019, Fig. C.1 -- Load Spectra and S-N Curve (Ref.)')
            
            legend({'SN', '\sigma_{HG}', '\sigma_{H1}'}, 'location', 'best');
            
            fig_axes = findobj(gcf, 'Type', 'Axes');
            set(fig_axes, 'fontName', 'Times');
            set(fig_axes, 'fontSize', 12);
            set(fig_axes, 'box', 'on');
            set(fig_axes, 'xScale', 'log');
            set(fig_axes, 'yScale', 'log');
            set(fig_axes, 'xlim', [1e2 1e10]);
            set(fig_axes, 'ylim', [800.0 3.0e3]);
            
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
        
        function val = dynamic_factor(obj)
            z1 = obj.z(1);
            uu = obj.u;
            cond = (obj.v_pitch_line*z1/100.0)*sqrt((uu.^2)/(1.0 + uu.^2));
            if(cond < 3.0) % [m/s]
                warning('MATISO_6336:KV', ['Calculating K_v using method B', ...
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
        
    end
    
    %% Set methods:
    methods
    end
    
    %% Get methods:
    methods
        function val = get.S_H(obj)
            val = obj.S_Hmin.*obj.sigma_HP./obj.sigma_H;
        end
        
        function val = get.S_F(obj)
            val = nan(size(obj.S_H));
        end
        
        function val = get.sigma_H0(obj)
            % Nominal contact stress at pitch point:
            num = obj.F_t*(obj.u + 1.0);
            den = obj.d(1)*obj.b*obj.u;
            val = obj.Z_H  * ...
                  obj.Z_E  * ...
                  obj.Z_eps* ...
                  obj.Z_beta*sqrt(num/den);
        end
        
        function val = get.sigma_H(obj)
            % Contact stress:
            val = obj.sigma_H0*obj.Z_BD*sqrt(obj.K_gamma* ...
                                             obj.K_A    * ...
                                             obj.K_v    * ...
                                             obj.K_Hbeta* ...
                                             obj.K_Halpha); % pinion/wheel
        end
        
        function val = get.sigma_HP(obj)
            % Permissible contact stress:
            val = obj.sigma_Hlim.*obj.Z_NT* ...
                                  obj.Z_L * ...
                                  obj.Z_v * ...
                                  obj.Z_R * ...
                                  obj.Z_W * ...
                                  obj.Z_X./obj.S_Hmin;
        end
        
        function val = get.sigma_Hlim(obj)
            val = [obj.material(1:2).sigma_Hlim];
        end
        
        function val = get.F_t(obj)
            % [N], Nominal tangential load:
            val = 2.0e3*(obj.T_nominal(1)/obj.d(1))/obj.N_p;
            val = abs(val);
        end
        
        function val = get.Z_NT(obj)
            val = zeros(1, 2);
            % assuming optimum lubrification, material, manufacturing and experience
%             y_lim = 1.0;
            y_lim = 0.85;
            
            for idx = 1:2
                line = obj.material(idx).row;

                switch line
                    case 1
                        % St, V, GGG (perl. bai.), GTS (perl.), Eh, IF (when limited pitting is permitted)
                        x = [6.0e5 1.0e7 1.0e9 1.0e10];
                        y = [1.6   1.3   1.0   y_lim]; 
                    case 2
                        % St, V, GGG (perl. bai.), GTS (perl.), Eh, IF
                        x = [1.0e5 5.0e7 1.0e10];
                        y = [1.6   1.0   y_lim];
                    case 3
                        % GG, GGG (ferr.), NT (nitr.), NV (nitr.)
                        x = [1.0e5 2.0e6 1.0e10];
                        y = [1.3   1.0   y_lim];
                    case 4
                        % NV (nitrocar.)
                        x = [1.0e5 2.0e6 1.0e10];
                        y = [1.1   1.0   y_lim];
                    otherwise
                        error("MATISO_6336:ZNT", "Invalid input [%d].\nValid options are 1 to 4.\n", line);
                end

                jdx = find(y == 1, 1, 'first');
                N_Lref = x(jdx);
                
                N = obj.N_L(idx);
                
                if(N <= x(1))
                    val(idx) = y(1);
                elseif(N > x(end))
                    val(idx) = y(end);
                else
                    if(N > N_Lref)
                        x = x([jdx end]);
                        y = y([jdx end]);
                    end
                    val(idx) = exp(interp1(log(x), log(y), log(N)));
                end
            end
            
        end
        
        function val = get.N_L(obj)
            % Number of load cycles:
            if(strcmp(obj.configuration, 'parallel'))
                N_L1 = obj.n_input*60.0*obj.L_h; % pinion
                val = N_L1*[1.0 1.0/obj.u];
            elseif(strcmp(obj.configuration, 'planetary'))
                n_rel = obj.n_nominal(1) - obj.n_nominal(4);
                val = [obj.N_p obj.z(1)/obj.z(2)]*60.0*obj.L_h*n_rel;
            end
        end
        
        function val = get.K_v(obj)
            val = obj.K_v_val;
        end
        
        function val = get.line_load(obj)
            val = obj.F_t*obj.K_A*obj.K_gamma/mean(obj.b);
        end
        
        function val = get.v_pitch_line(obj)
            val = nan;
            if(strcmp(obj.configuration, 'parallel'))
                val = (pi*obj.n_input/60.0e3)*obj.d(1);
            elseif(strcmp(obj.configuration, 'planetary'))
                val = (pi*obj.n_output/60.0e3)*(obj.d(1) - obj.a_w/obj.u);
            end
        end
        
        function val = get.K_Halpha(obj)
            val = obj.K_Falpha;
            
            % Section 8.3.4 of [1]:
            % Transverse load factor (contact stress), Eq. (75):
            K_Halpha_lim = obj.eps_gamma/(obj.eps_alpha*obj.Z_eps.^2);
            if(val > K_Halpha_lim)
                val = K_Halpha_lim;
            elseif(val < 1.0)
                val = 1.0;
            end            
        end
        
        function val = get.K_Falpha(obj)
            term = obj.c_gamma_alpha*(obj.f_pb - obj.y_alpha)/(obj.F_tH/obj.b);
            % Section 8.3.2 of [1]:
            if(obj.eps_gamma <= 2.0) % Eq. (73)
                val = (0.9 + 0.4*term)*(obj.eps_gamma/2.0);
            else % Eq. (74)
                val = 0.9 + 0.4*sqrt(2.0*(obj.eps_gamma - 1.0)/obj.eps_gamma)*term;
            end
            
            % Section 8.3.5 of [1]:
            % Transverse load factor (root stress), Eq. (76):
            K_Falpha_lim = obj.eps_gamma/(0.25*obj.eps_alpha + 0.75);
            if(val > K_Falpha_lim)
                val = K_Falpha_lim;
            elseif(val < 1.0)
                val = 1.0;
            end
        end
        
        function val = get.K_Hbeta(obj)
            % Alternative: according to [5], Eq.(17.14):
            % val = 1.0 + 0.4*(obj.b/obj.d(1))^2;
            
            % Mean transverse tangential load at the reference circle:
            F_m = obj.F_t*obj.K_A*obj.K_gamma*obj.K_v;
            
            k_unit = 0.023; % [mm-um/N]
            B_star = 1.0 + 2.0*(obj.N_p - 1.0);
            
            % Component of equivalent misalignment due to deformations of
            % pinion and wheel shafts. 
            % Assumption made: Pinion offset is null, s = 0.
            f_sh = (F_m/mean(obj.b))*k_unit*(abs(B_star - 0.3) + 0.3)*(mean(obj.b)/obj.d(1))^2;
            
            % mesh misalignment due to manufacturing deviations:
            f_ma = norm(obj.f_Hbeta);
            
            % Initial equivalent misalignment (before running-in), assuming
            % that no helix modification was applied:
            F_betax = 1.33*f_sh + f_ma;
            
            % Running-in allowance:
            y_beta = 0.15*F_betax;
            
            F_betay = F_betax - y_beta;
            
            cond = (F_betay*obj.c_gamma_beta)/(2.0*F_m/mean(obj.b));
            
            if(cond >= 1.0)
                val = 2.0*sqrt(cond);
            else
                val = 1.0 + cond;
            end
            
        end
        
        function val = get.K_Fbeta(obj)
            % Method B, Sec. 7.6:
            % h_1 = h_aP + h_fP + k_1*m_n;
            h_1 = abs(obj.d_a(1) - obj.d_f(1))/2.0;
            bh1 = obj.b/h_1;
            h_2 = abs(obj.d_a(2) - obj.d_f(2))/2.0;
            bh2 = obj.b/h_2;
            
            bh = min(bh1, bh2);
            
            if(bh < 3.0)
                bh = 3.0;
            end
            
            N_F = (bh.^2)/(1.0 + bh + bh.^2);
            
            val = power(obj.K_Hbeta, N_F);
        end
        
        function val = get.Z_E(obj)
            % [N/mm^2], Young's modulus:
            E = [obj.material.E];
            % [-], Poisson's ratio:
            nu = [obj.material.nu];
            
            vec = (1.0 - nu.^2)./E;
            val = sqrt(1.0/(pi*sum(vec(1:2))));
            
        end
        
        function val = get.Z_H(obj)
            num = 2.0*cosd(obj.beta_b)*cosd(obj.alpha_wt);
            den = sind(obj.alpha_wt)*cosd(obj.alpha_t).^2;
            val = sqrt(num/den);
        end
        
        function val = get.Z_eps(obj)
            eps_a = obj.eps_alpha;
            eps_b = obj.eps_beta;
            
            if(obj.beta == 0.0)
                val = sqrt((4.0 - eps_a)/3.0);
            else
                if(eps_b < 1.0)
                    val = sqrt((1.0 - eps_b)*(4.0 - eps_a)/3.0 + eps_b/eps_a);
                else
                    val = sqrt(1.0/eps_a);
                end
            end
            
        end
        
        function val = get.Z_beta(obj)
            % Helix angle factor: (sec. 9)
            val = 1.0/sqrt(cosd(obj.beta));
        end
        
        function val = get.Z_BD(obj)
            M_1 = tand(obj.alpha_wt)/sqrt((sqrt((obj.d_a(1)/obj.d_b(1))^2 - 1.0) - ...
                2.0*pi/obj.z(1))*(sqrt((obj.d_a(2)/obj.d_b(2))^2 - 1.0) - ...
                (obj.eps_alpha - 1.0)*2.0*pi/obj.z(2)));
            M_2 = tand(obj.alpha_wt)/sqrt((sqrt((obj.d_a(2)/obj.d_b(2))^2 - 1.0) - ...
                2.0*pi/obj.z(2))*(sqrt((obj.d_a(1)/obj.d_b(1))^2 - 1.0) - ...
                (obj.eps_alpha - 1.0)*2.0*pi/obj.z(1)));
            M = [M_1 M_2];
            
            val = nan(1, 2);
            if((obj.eps_beta == 0.0) && (obj.eps_alpha > 1.0))
                for idx = 1:2
                    if(M(idx) > 1.0)
                        val(idx) = M(idx);
                    else
                        val(idx) = 1.0;
                    end
                end
            elseif((obj.eps_alpha > 1.0) && (obj.eps_beta >= 1.0))
                val = ones(1, 2);
            elseif((obj.eps_alpha > 1.0) && (obj.eps_beta <  1.0))
                val = M - obj.eps_beta*(M - 1.0);
            end
        end
        
        function val = get.Z_R(obj)
            rho_1 = 0.5*obj.d_b(1)*tand(obj.alpha_wt);
            rho_2 = 0.5*obj.d_b(2)*tand(obj.alpha_wt);
            
            rho_red = (rho_1*rho_2)/(rho_1 + rho_2);
            
            R_z10 = obj.R_z*power(10.0/rho_red, 1.0/3.0);
            sig_Hlim = min(obj.sigma_Hlim);
            
            if(sig_Hlim < 850.0) % [N/mm^2]
                C_ZR = 0.15;
            elseif((850.0 <= sig_Hlim) && (sig_Hlim < 1200.0))
                C_ZR = 0.32 - sig_Hlim*2.0e-4;
            else
                C_ZR = 0.08;
            end
            
            val = power(3.0/R_z10, C_ZR);
        end
        
        function val = get.Z_L(obj)
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
            val = C_ZL + 4.0*(1.0 - C_ZL)/(1.2 + 134.0/obj.nu_40)^2;
            
        end
        
        function val = get.Z_v(obj)
            % [N/mm^2],  Allowable contact stress number:
            sig_Hlim = min(obj.sigma_Hlim);
            
            if(sig_Hlim  < 850.0) % [N/mm^2]
                C_ZL = 0.83;
            elseif((850.0 <= sig_Hlim) && (sig_Hlim  < 1200.0))
                C_ZL = sig_Hlim/4375.0 + 0.6357;
            else
                C_ZL = 0.91;
            end
            
            % Velocity factor:
            C_Zv = C_ZL + 0.02;
            val = C_Zv + 2.0*(1.0 - C_Zv)/sqrt(0.8 + 32.0/obj.v_pitch_line);
            
        end
        
        function val = get.Z_W(~)
            val = 1.0;
        end
        
        function val = get.Z_X(~)
            val = 1.0;
        end
        
    end
end
