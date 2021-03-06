classdef Lin_Parker_99 < Dynamic_Formulation
    %LIN_PARKER_99 Calculates the inertia and stiffness matrices of a 
    % multi-stage Drivetrain object according to:
    % [1] J. Lin and R. G. Parker, 'Analytical Characterization of the
    % Unique Properties of Planetary Gear Free Vibration', Journal
    % of Vibration and Acoustics, vol. 121, no. 3, pp. 316-321,
    % 1999. https://doi.org/10.1115/1.2893982
    % [2] J. Lin, 'Analytical investigation of planetary gear dynamics',
    % Ph.D thesis, Ohio State University, 2000.
    % http://rave.ohiolink.edu/etdc/view?acc_num=osu1488203552779634
    % [3] C. G. Cooley and R. G. Parker, 'Vibration Properties of
    % High-Speed Planetary Gears With Gyroscopic Effects' Journal of
    % Vibration and Acoustics, vol. 134, no. 6, Dec. 2012.
    % https://doi.org/10.1115/1.4006646
    %
    
    properties
        D_bearing; % Bearing damping
        D_mesh;    % Mesh damping
        
    end
    
    properties(Dependent)
        KK;
    end
    
    methods
        function obj = Lin_Parker_99(DT, varargin)
            default = {'fault_type',  '', ...
                       'fault_val' , 0.0};
            
            default = scaling_factor.process_varargin(default, varargin);
            
            obj@Dynamic_Formulation(DT, 'fault_type', default.fault_type, ...
                                        'fault_val' , default.fault_val);
        end
        
        %% Calculation:
        function [MM, GG] = inertia_matrix(obj)
            N = obj.n_DOF(end);
            MM = zeros(N, N);
            GG = zeros(N, N);
            
            DT = obj.drive_train;

            m_r = DT.m_Rotor;            m_g = DT.m_Gen;
            J_r = DT.J_Rotor;            J_g = DT.J_Gen;
            
            range = 1:3;
            MM(range, range) = MM(range, range) + diag([m_r, m_r, J_r]);
            
            range = (N-2):N;
            MM(range, range) = MM(range, range) + diag([m_g, m_g, J_g]);
            
            M_tmp = DT.main_shaft.inertia_matrix('Lin_Parker_99');
            
            range = 1:6;
            MM(range, range) = MM(range, range) + M_tmp;
            
            N = obj.n_DOF;
            Z3 = zeros(3);
            for idx = 1:DT.N_stage
                R = Lin_Parker_99.stage_coordinate_change(DT.stage(idx));
                [M_idx, G_idx] = Lin_Parker_99.stage_inertia_matrix(DT.stage(idx));
                M_idx = R' * M_idx * R;
                G_idx = R' * G_idx * R;
                
                range = (N(idx) - 2):(N(idx + 1));
                
                % adding shaft inertia:
                n_idx = length(range);
                rng = (-5:0) + n_idx;
                
                M_tmp = blkdiag(M_idx, Z3);
                M_tmp(rng, rng) = M_tmp(rng, rng) + ...
                    DT.stage(idx).output_shaft.inertia_matrix('Lin_Parker_99');
                
                MM(range, range) = MM(range, range) + M_tmp;
                GG(range, range) = GG(range, range) + blkdiag(G_idx, Z3);
            end
        end
        
        function D = damping_matrix(obj)
            N = obj.n_DOF(end);
            Db = zeros(N, N);
            Dm = zeros(N, N);
            
            DT = obj.drive_train;
            
            D_tmp = DT.main_shaft.damping_matrix('Lin_Parker_99');
            
            range = 1:6;
            Db(range, range) = Db(range, range) + D_tmp;
            
            N = obj.n_DOF;
            Z3 = zeros(3);
            for idx = 1:DT.N_stage
                R = Lin_Parker_99.stage_coordinate_change(DT.stage(idx));
                [Db_idx, Dm_idx] = Lin_Parker_99.stage_stiffness_matrix(DT.stage(idx));
                Db_idx = R' * Db_idx * R;
                Dm_idx = R' * Dm_idx * R;
                
                range = (N(idx) - 2):(N(idx + 1));
                
                % adding shaft stiffness:
                n_idx = length(range);
                rng = (-5:0) + n_idx;
                
                Db_tmp = blkdiag(Db_idx, Z3);
                Db_tmp(rng, rng) = Db_tmp(rng, rng) + ...
                    DT.stage(idx).output_shaft.damping_matrix('Lin_Parker_99');
                
                Db(range, range) = Db(range, range) + Db_tmp;
                Dm(range, range) = Dm(range, range) + blkdiag(Dm_idx, Z3);
            end
            D = Dm + Db;
        end
        
        function [Km, Kb, KOm] = stiffness_matrix(obj)
            N = obj.n_DOF(end);
            Kb  = zeros(N);     Km  = zeros(N);
            KOm = zeros(N);
            
            DT = obj.drive_train;

            K_shaft = DT.main_shaft.stiffness_matrix('Lin_Parker_99');
            
            range = 1:6;
            Kb(range, range) = Kb(range, range) + K_shaft;

            N = obj.n_DOF;
            Z3 = zeros(3);
            for idx = 1:DT.N_stage
                [Km_idx, Kb_idx, KO_idx] = Lin_Parker_99.stage_stiffness_matrix(DT.stage(idx));
                
                R = Lin_Parker_99.stage_coordinate_change(DT.stage(idx));
                Km_idx = R' * Km_idx * R;
                Kb_idx = R' * Kb_idx * R;
                KO_idx = R' * KO_idx * R;
                
                range = (N(idx) - 2):(N(idx + 1));

                % adding shaft stiffness:
                n_idx = length(range);
                rng = (-5:0) + n_idx;
                
                Kb_tmp = blkdiag(Kb_idx, Z3);
                Kb_tmp(rng, rng) = Kb_tmp(rng, rng) + ...
                    DT.stage(idx).output_shaft.stiffness_matrix('Lin_Parker_99');

                Km(range, range)  = Km( range, range) + blkdiag(Km_idx, Z3);
                Kb(range, range)  = Kb( range, range) + Kb_tmp;
                KOm(range, range) = KOm(range, range) + blkdiag(KO_idx, Z3);
            end
        end
        
        function cc = centripetal_force_vector(obj)
            N = obj.n_DOF;
            cc = zeros(N(end), 1);
            DT = obj.drive_train;
            for idx = 1:DT.N_stage
                R = Lin_Parker_99.stage_coordinate_change(DT.stage(idx));
                cc_idx = Lin_Parker_99.stage_centripetal_force_vector(DT.stage(idx));
                cc_idx = R' * cc_idx;
                
                range = (N(idx) - 2):(N(idx + 1) - 3);
                cc(range) = cc(range) + cc_idx;
            end
        end
        
        function bb = external_load_vector(obj)
            bb         = zeros(obj.n_DOF(end), 2);
            bb(3  , 1) = 1.0; % Torque at the rotor 
            bb(end, 2) = 1.0; % Torque at the generator
        end
        
        function nn = calculate_DOF(obj)
            DT = obj.drive_train;

            nn = ones(DT.N_stage + 1, 1)*6.0;
            
            for idx = 1:DT.N_stage
                if(strcmp(DT.stage(idx).configuration, 'parallel'))
                    tmp = DT.stage(idx).N_p + 1;
                elseif(strcmp(DT.stage(idx).configuration, 'planetary'))
                    tmp = DT.stage(idx).N_p + 2;
                end
                
                nn(idx + 1) = nn(idx) + tmp*3;
            end
        end
        
    end
    
    methods(Static)
        function [MM, GG] = stage_inertia_matrix(stage_idx)
            G_ = @(x)[0.0  , -2.0*x, 0.0;
                      2.0*x,  0.0  , 0.0;
                      0.0  ,  0.0  , 0.0];
            
            if(strcmp(stage_idx.configuration, 'parallel'))
                m_p = stage_idx.mass(1);
                m_w = stage_idx.mass(2);
                
                r_p = (stage_idx.d_b(1)*1.0e-3)/2.0;
                r_w = (stage_idx.d_b(2)*1.0e-3)/2.0;
                
                J_p = stage_idx.J_x(1)/r_p^2;
                J_w = stage_idx.J_x(2)/r_w^2;
                
                diag_inertia = [m_p, m_p, J_p, ... % pinion
                                m_w, m_w, J_w];    % wheel
                
                % Gyroscopic matrix:
                GG = blkdiag(G_(m_p), G_(m_w));
                
            elseif(strcmp(stage_idx.configuration, 'planetary'))
                % carrier, ring, sun, N planets
                
                m_c = stage_idx.carrier.mass;
                m_r = stage_idx.mass(3);
                m_s = stage_idx.mass(1);
                m_p = stage_idx.mass(2);
                
                r_c =  stage_idx.a_w   *1.0e-3;
                r_r = (stage_idx.d_b(3)*1.0e-3)/2.0;
                r_s = (stage_idx.d_b(1)*1.0e-3)/2.0;
                r_p = (stage_idx.d_b(2)*1.0e-3)/2.0;
                
                J_c = stage_idx.carrier.J_x/r_c^2;
                J_r = stage_idx.J_x(3)     /r_r^2;
                J_s = stage_idx.J_x(1)     /r_s^2;
                J_p = stage_idx.J_x(2)     /r_p^2;
                
                diag_inertia = [m_c, m_c, J_c, ...
                                m_r, m_r, J_r, ...
                                m_s, m_s, J_s, ...
                               [m_p, m_p, J_p]*repmat(eye(3), 1, stage_idx.N_p)];
                
                % Gyroscopic matrix:
                Gp = repmat({G_(m_p)}, 1, stage_idx.N_p);
                Gp = blkdiag(Gp{:});
                
                GG = blkdiag(G_(m_c), ... % carrier
                             G_(m_r), ... % ring
                             G_(m_s), ... % sun
                             Gp     );    % planet
                
            end
            
            MM = diag(diag_inertia);
        end
        
        function [Db, Dm] = stage_damping_matrix(stage_idx)
            alpha_n = stage_idx.alpha_n;
            
            N1    = stage_idx.N_p - 1;
            psi   = (0:N1)*(360.0/stage_idx.N_p);
            psi_s = psi - alpha_n;
            psi_r = psi + alpha_n;
            
            % carrier-carrier bearing damping sub-matrix:
            D_c1 = @(i)[ 1.0         , 0.0         , -sind(psi(i));
                         0.0         , 1.0         ,  cosd(psi(i));
                        -sind(psi(i)), cosd(psi(i)),        1.0  ];

            % carrier-planet bearing damping sub-matrix:
            D_c2 = @(i)[-cosd(psi(i)),  sind(psi(i)), 0.0;
                        -sind(psi(i)), -cosd(psi(i)), 0.0;
                         0.0         , -1.0         , 0.0];
            
            % ring-ring mesh damping sub-matrix:
            D_r1 = @(i)[ sind(psi_r(i))^2             , -sind(psi_r(i))*cosd(psi_r(i))  , -sind(psi_r(i));
                        -sind(psi_r(i))*cosd(psi_r(i)),                 cosd(psi_r(i))^2,  cosd(psi_r(i));
                        -sind(psi_r(i))               ,                 cosd(psi_r(i))  ,       1.0     ];
            
            % ring-planet mesh damping sub-matrix:
            D_r2 = @(i)[-sind(psi_r(i))*sind(alpha_n),  sind(psi_r(i))*cosd(alpha_n),  sind(psi_r(i));
                         cosd(psi_r(i))*sind(alpha_n), -cosd(psi_r(i))*cosd(alpha_n), -cosd(psi_r(i));
                                        sind(alpha_n), -               cosd(alpha_n), -1.0          ];
            
            % planet-ring mesh damping sub-matrix:
            D_r3 = [ sind(alpha_n)^2            , -sind(alpha_n)*cosd(alpha_n)  , -sind(alpha_n);
                    -sind(alpha_n)*cosd(alpha_n),                cosd(alpha_n)^2,  cosd(alpha_n);
                    -sind(alpha_n)              ,                cosd(alpha_n)  ,        1.0   ];
            
            % sun-sun mesh-damping matrix:
            D_s1 = @(i)[                sind(psi_s(i))^2, -cosd(psi_s(i))*sind(psi_s(i)), -sind(psi_s(i));
                        -cosd(psi_s(i))*sind(psi_s(i))  ,  cosd(psi_s(i))^2             ,  cosd(psi_s(i));
                        -               sind(psi_s(i))  ,  cosd(psi_s(i))               ,         1.0   ];
            
            % sun-planet mesh-damping matrix:
            D_s2 = @(i)[ sind(psi_s(i))*sind(alpha_n),  sind(psi_s(i))*cosd(alpha_n), -sind(psi_s(i));
                        -cosd(psi_s(i))*sind(alpha_n), -cosd(psi_s(i))*cosd(alpha_n),  cosd(psi_s(i));
                        -               sind(alpha_n), -               cosd(alpha_n),         1.0   ];
                                 
            % planet-planet (?) mesh-damping matrix:
            D_s3 = [ sind(alpha_n)^2            ,  sind(alpha_n)*cosd(alpha_n)  , -sind(alpha_n);
                     sind(alpha_n)*cosd(alpha_n),                cosd(alpha_n)^2, -cosd(alpha_n);
                    -sind(alpha_n)              , -              cosd(alpha_n)  ,        1.0   ];
            
            if(strcmp(stage_idx.configuration, 'parallel'))
                
                r_p = (stage_idx.d_b(1)*1.0e-3)/2.0; % [m]
                r_w = (stage_idx.d_b(2)*1.0e-3)/2.0;
                
                % Bearing component:
                brg_w  = stage_idx.bearing(4:6);
                brg_w = brg_w.parallel_association();
                d_wx = brg_w.C_y;
                d_wy = brg_w.C_z;
                d_wu = brg_w.C_alpha/(r_w^2);
                
                Db = diag([zeros(1, 3) d_wx d_wy d_wu]);
                
                brg_p  = stage_idx.bearing(1:3);
                brg_p = brg_p.parallel_association();
                d_px = brg_p.C_y;
                d_py = brg_p.C_z;
                d_pu = brg_p.C_alpha/(r_p^2);
                
                % Mesh component:
                d_pw = 500.0e6;
                D_c31 =  diag([d_px, d_py, d_pu]);
                Dm = [d_pw*D_s3 + D_c31, d_pw*D_s2(1) ;
                      d_pw*D_s2(1)'    , d_pw*D_s1(1)];

            elseif(strcmp(stage_idx.configuration, 'planetary'))
                np = 3*stage_idx.N_p;
                
                r_c =  stage_idx.a_w   *1.0e-3; % [m]
                r_r = (stage_idx.d_b(3)*1.0e-3)/2.0;
                r_s = (stage_idx.d_b(1)*1.0e-3)/2.0;
                r_p = (stage_idx.d_b(2)*1.0e-3)/2.0;
                
                % Bearing component:
                brg_s = stage_idx.bearing(5);
                d_sx = brg_s.C_y;
                d_sy = brg_s.C_z;
                d_su = brg_s.C_alpha/(r_s^2);

                brg_r = stage_idx.bearing(6);
                d_rx = brg_r.C_y;
                d_ry = brg_r.C_z;
                d_ru = brg_r.C_alpha/(r_r^2);
                
                brg_c = stage_idx.bearing(3:4);
                brg_c = brg_c.parallel_association();
                d_cx = brg_c.C_y;
                d_cy = brg_c.C_z;
                d_cu = brg_c.C_alpha/(r_c^2);
                
                Db = diag([d_cx, d_cy, d_cu, ... % carrier
                           d_rx, d_ry, d_ru, ... % ring
                           d_sx, d_sy, d_su, ... % sun
                           zeros(1, np)]);       % planets
                
                brg_p = stage_idx.bearing(1:2);
                brg_p = brg_p.parallel_association();
                d_px = brg_p.C_y;
                d_py = brg_p.C_z;
                d_pu = brg_p.C_alpha/(r_p^2);

                d_sp = 500.0e6;
                d_pr = 500.0e6;
                
                D_c_row = zeros(3, np);
                D_r_row = D_c_row;
                D_s_row = D_c_row;
                
                sum_Dc1 = zeros(3);
                sum_Dr1 = zeros(3);
                sum_Ds1 = zeros(3);
                for idx = 1:stage_idx.N_p
                    rng = (-2:0) + 3*idx;
                    D_c_row(:, rng) = D_c2(idx);
                    D_r_row(:, rng) = D_r2(idx);
                    D_s_row(:, rng) = D_s2(idx);
                    
                    sum_Dc1 = sum_Dc1 + D_c1(idx);
                    sum_Dr1 = sum_Dr1 + D_r1(idx);
                    sum_Ds1 = sum_Ds1 + D_s1(idx);
                end
                
                D_c_row = diag([d_px, d_py, d_pu])*D_c_row;
                D_r_row =                d_pr  *D_r_row;
                D_s_row =                d_sp  *D_s_row;
                sum_Dc1 = diag([d_px, d_py, d_pu])*sum_Dc1;
                sum_Dr1 =                d_pr  *sum_Dr1;
                sum_Ds1 =                d_sp  *sum_Ds1;
                
                diag_01 = blkdiag(sum_Dc1, sum_Dr1, sum_Ds1);
                diag_up = [D_c_row;
                           D_r_row;
                           D_s_row];

                D_pp = diag([d_px, d_py, d_pu]) + d_pr*D_r3 + d_sp*D_s3;
                DppCell = repmat({D_pp}, 1, stage_idx.N_p);
                diag_02 = blkdiag(DppCell{:});
                
                Dm = [diag_01 , diag_up;
                      diag_up', diag_02];
            end
            
            % removing rounding errors from matrices:
            Db(abs(Db) < 1.0e-4) = 0.0;
            Dm(abs(Dm) < 1.0e-4) = 0.0;
            
        end
        
        function [Km, Kb, KOm] = stage_stiffness_matrix(stage_idx)
            alpha_n = stage_idx.alpha_n;
            
            N1    = stage_idx.N_p - 1;
            psi   = (0:N1)*(360.0/stage_idx.N_p);
            psi_s = psi - alpha_n;
            psi_r = psi + alpha_n;
            
            % carrier-carrier bearing stiffness sub-matrix:
            K_c1 = @(i)[ 1.0         , 0.0         , -sind(psi(i));
                         0.0         , 1.0         ,  cosd(psi(i));
                        -sind(psi(i)), cosd(psi(i)),        1.0  ];

            % carrier-planet bearing stiffness sub-matrix:
            K_c2 = @(i)[-cosd(psi(i)),  sind(psi(i)), 0.0;
                        -sind(psi(i)), -cosd(psi(i)), 0.0;
                         0.0         , -1.0         , 0.0];
            
            % ring-ring mesh stiffness sub-matrix:
            K_r1 = @(i)[ sind(psi_r(i))^2             , -sind(psi_r(i))*cosd(psi_r(i))  , -sind(psi_r(i));
                        -sind(psi_r(i))*cosd(psi_r(i)),                 cosd(psi_r(i))^2,  cosd(psi_r(i));
                        -sind(psi_r(i))               ,                 cosd(psi_r(i))  ,       1.0     ];
            
            % ring-planet mesh stiffness sub-matrix:
            K_r2 = @(i)[-sind(psi_r(i))*sind(alpha_n),  sind(psi_r(i))*cosd(alpha_n),  sind(psi_r(i));
                         cosd(psi_r(i))*sind(alpha_n), -cosd(psi_r(i))*cosd(alpha_n), -cosd(psi_r(i));
                                        sind(alpha_n), -               cosd(alpha_n), -1.0          ];
            
            % planet-ring mesh stiffness sub-matrix:
            K_r3 = [ sind(alpha_n)^2            , -sind(alpha_n)*cosd(alpha_n)  , -sind(alpha_n);
                    -sind(alpha_n)*cosd(alpha_n),                cosd(alpha_n)^2,  cosd(alpha_n);
                    -sind(alpha_n)              ,                cosd(alpha_n)  ,        1.0   ];
            
            % sun-sun mesh-stiffness matrix:
            K_s1 = @(i)[                sind(psi_s(i))^2, -cosd(psi_s(i))*sind(psi_s(i)), -sind(psi_s(i));
                        -cosd(psi_s(i))*sind(psi_s(i))  ,  cosd(psi_s(i))^2             ,  cosd(psi_s(i));
                        -               sind(psi_s(i))  ,  cosd(psi_s(i))               ,         1.0   ];
            
            % sun-planet mesh-stiffness matrix:
            K_s2 = @(i)[ sind(psi_s(i))*sind(alpha_n),  sind(psi_s(i))*cosd(alpha_n), -sind(psi_s(i));
                        -cosd(psi_s(i))*sind(alpha_n), -cosd(psi_s(i))*cosd(alpha_n),  cosd(psi_s(i));
                        -               sind(alpha_n), -               cosd(alpha_n),         1.0   ];
                                 
            % planet-planet (?) mesh-stiffness matrix:
            K_s3 = [ sind(alpha_n)^2            ,  sind(alpha_n)*cosd(alpha_n)  , -sind(alpha_n);
                     sind(alpha_n)*cosd(alpha_n),                cosd(alpha_n)^2, -cosd(alpha_n);
                    -sind(alpha_n)              , -              cosd(alpha_n)  ,        1.0   ];
            
            if(strcmp(stage_idx.configuration, 'parallel'))
                % Centripetal component:
                m_1 = stage_idx.mass(1);
                m_2 = stage_idx.mass(2);
                KOm = diag([m_1, m_1, 0, ... % pinion
                            m_2, m_2, 0]);   % wheel
                
                r_p = (stage_idx.d_b(1)*1.0e-3)/2.0; % [m]
                r_w = (stage_idx.d_b(2)*1.0e-3)/2.0;
                
                % Bearing component:
                brg_w  = stage_idx.bearing(4:6);
                brg_w = brg_w.parallel_association();
                k_wx = brg_w.K_y;
                k_wy = brg_w.K_z;
                k_wu = brg_w.K_alpha/(r_w^2);
                
                Kb = diag([zeros(1, 3) k_wx k_wy k_wu]);
                
                brg_p  = stage_idx.bearing(1:3);
                brg_p = brg_p.parallel_association();
                k_px = brg_p.K_y;
                k_py = brg_p.K_z;
                k_pu = brg_p.K_alpha/(r_p^2);
                
                % Mesh component:
                k = stage_idx.k_mesh;
                K_c31 =  diag([k_px, k_py, k_pu]);
                Km = [k*K_s3 + K_c31, k*K_s2(1) ;
                      k*K_s2(1)'    , k*K_s1(1)];

            elseif(strcmp(stage_idx.configuration, 'planetary'))
                np = 3*stage_idx.N_p;
                
                % Centripetal component:
                m_c = stage_idx.carrier.mass;
                m_r = stage_idx.mass(3);
                m_s = stage_idx.mass(1);
                m_p = stage_idx.mass(2);
                                                          
                KOm = diag([m_c m_c 0.0, ... % carrier
                            m_r m_r 0.0, ... % ring
                            m_s m_s 0.0, ... % sun
                           [m_p m_p 0.0]*repmat(eye(3), 1, stage_idx.N_p)]);
                
                r_c =  stage_idx.a_w   *1.0e-3; % [m]
                r_r = (stage_idx.d_b(3)*1.0e-3)/2.0;
                r_s = (stage_idx.d_b(1)*1.0e-3)/2.0;
                r_p = (stage_idx.d_b(2)*1.0e-3)/2.0;
                
                % Bearing component:
                brg_s = stage_idx.bearing(5);
                k_sx = brg_s.K_y;
                k_sy = brg_s.K_z;
                k_su = brg_s.K_alpha/(r_s^2);

                brg_r = stage_idx.bearing(6);
                k_rx = brg_r.K_y;
                k_ry = brg_r.K_z;
                k_ru = brg_r.K_alpha/(r_r^2);
                
                brg_c = stage_idx.bearing(3:4);
                brg_c = brg_c.parallel_association();
                k_cx = brg_c.K_y;
                k_cy = brg_c.K_z;
                k_cu = brg_c.K_alpha/(r_c^2);
                
                Kb = diag([k_cx, k_cy, k_cu, ... % carrier
                           k_rx, k_ry, k_ru, ... % ring
                           k_sx, k_sy, k_su, ... % sun
                           zeros(1, np)]);       % planets
                
                brg_p = stage_idx.bearing(1:2);
                brg_p = brg_p.parallel_association();
                k_px = brg_p.K_y;
                k_py = brg_p.K_z;
                k_pu = brg_p.K_alpha/(r_p^2);

                % Mesh component:
                sun_pla = stage_idx.sub_set('sun_planet');
                pla_rng = stage_idx.sub_set('planet_ring');
                
                % OBS.: got a negative k_mesh for DTU_10MW planet-ring
                k_sp = sun_pla.k_mesh;
                k_pr = pla_rng.k_mesh;
                
                K_c_row = zeros(3, np);
                K_r_row = K_c_row;
                K_s_row = K_c_row;
                
                sum_Kc1 = eye(3)*stage_idx.N_p;
                sum_Kr1 = diag([1.0/2.0, 1.0/2.0, 1.0])*stage_idx.N_p;
                sum_Ks1 = diag([1.0/2.0, 1.0/2.0, 1.0])*stage_idx.N_p;
                for idx = 1:stage_idx.N_p
                    rng = (-2:0) + 3*idx;
                    K_c_row(:, rng) = K_c2(idx);
                    K_r_row(:, rng) = K_r2(idx);
                    K_s_row(:, rng) = K_s2(idx);
                end
                
                K_c_row = diag([k_px, k_py, k_pu])*K_c_row;
                K_r_row =                   k_pr  *K_r_row;
                K_s_row =                   k_sp  *K_s_row;
                sum_Kc1 = diag([k_px, k_py, k_pu])*sum_Kc1;
                sum_Kr1 =                   k_pr  *sum_Kr1;
                sum_Ks1 =                   k_sp  *sum_Ks1;
                
                diag_01 = blkdiag(sum_Kc1, sum_Kr1, sum_Ks1);
                diag_up = [K_c_row;
                           K_r_row;
                           K_s_row];

                K_pp = diag([k_px, k_py, 0.0]) + k_pr*K_r3 + k_sp*K_s3;
                KppCell = repmat({K_pp}, 1, stage_idx.N_p);
                diag_02 = blkdiag(KppCell{:});
                
                Km = [diag_01 , diag_up;
                      diag_up', diag_02];
            end
            
            % removing rounding errors from matrices:
            Km( abs(Km)  < 1.0e-4) = 0.0;
            Kb( abs(Kb)  < 1.0e-4) = 0.0;
            KOm(abs(KOm) < 1.0e-4) = 0.0;
            
        end
        
        function cc = stage_centripetal_force_vector(stage_idx)
            if(strcmp(stage_idx.configuration, 'parallel'))
                cc = zeros(6, 1);
            elseif(strcmp(stage_idx.configuration, 'planetary'))
                n = 3.0*(stage_idx.N_p + 3.0);
                cc = zeros(n, 1);
                
                m_p = stage_idx.mass(2);
                r_s = (stage_idx.d_b(1)*1.0e-3)/2.0;
                r_p = (stage_idx.d_b(2)*1.0e-3)/2.0;
                
                cc(10:3:end) = m_p*(r_s + r_p);
            end
            
        end
        
        function [WVR, R, V, W] = stage_coordinate_change(stage_idx)
            %STAGE_COORDINATE_CHANGE performs the following coordinate
            %transformations:
            % R: to go from translational (x, y, u) to rotational
            % coordinates (x, y, theta), with u = r * theta.
            % V: to change element order, from:
            % - (Carrier, Sun, Planet) to (Carrier, Planet, Sun), or;
            % - (Pinion, Wheel) to (Wheel, Pinion)
            % W: to eliminate the ring gear DOFs.
            %
            
            Z3 = zeros(3);
            I3 = eye(3);
            
            if(strcmp(stage_idx.configuration, 'parallel'))
                r_p = (stage_idx.d_b(1)*1.0e-3)/2.0; % [m]
                r_w = (stage_idx.d_b(2)*1.0e-3)/2.0;
                
                diag_R = [1.0, 1.0, r_p, ...
                          1.0, 1.0, r_w];

                V = [Z3, I3;
                     I3, Z3];

                W = eye(6);
            elseif(strcmp(stage_idx.configuration, 'planetary'))
                np = 3*stage_idx.N_p;
                
                r_c =  stage_idx.a_w   *1.0e-3;      % [m]
                r_s = (stage_idx.d_b(1)*1.0e-3)/2.0;
                r_p = (stage_idx.d_b(2)*1.0e-3)/2.0;
                
                diag_R = [1.0, 1.0, r_c, ... % carrier
                         [1.0, 1.0, r_p]*repmat(eye(3), 1, stage_idx.N_p), ... % planet
                          1.0, 1.0, r_s]; % sun
                
                n = length(diag_R);
                W = eye(n);
                W = [W(1:3, :);
                     zeros(3, n);
                     W(4:end, :)];
                
                V = zeros(n);
                V(1:3, 1:3) = I3;
                
                col = np + 4;
                V(4:6, col:end) = I3;
                
                col = np + 3;
                V(7:end, 4:col) = eye(np);
            end
            
            R = diag(diag_R);
            WVR = W*V*R;
        end
        
        function flag = test_01(num_planet)
            % Testing the code using the parameters from:
            % [1], Table 1 and
            % [2], Table C.1;
            % Reference values for resonances were taken from:
            % [1], Table 2 and;
            % [2], Table 3.1;
            %
            % +-----------------------------+-------+-------+---------+--------+
            % |                             |  Sun  |  Ring | Carrier | Planet |
            % +-----------------------------+-------+-------+---------+--------+
            % |         Mass (kg)           |  0.4  |  2.35 |   5.43  |  0.66  |
            % +-----------------------------+-------+-------+---------+--------+
            % |         I/r^2 (kg)          |  0.39 |  3.00 |   6.29  |  0.61  |
            % +-----------------------------+-------+-------+---------+--------+
            % |     Base diameter (mm)      | 77.4  | 275.0 |  176.8  |  100.3 |
            % +-----------------------------+-------+-------+---------+--------+
            % |    Mesh stiffness (N/m)     |               5.0e8              |
            % +-----------------------------+-------+-------+---------+--------+
            % |  Bearing stiffness (N/m)    | 1.0e8 | 1.0e8 |  1.0e8  |  1.0e8 |
            % +-----------------------------+-------+-------+---------+--------+
            % |  Torsional stiffness (N/m)  |  0.0  | 1.0e9 |   0.0   |  1.0e8 |
            % +-----------------------------+-------+-------+---------+--------+
            % |    Pressure angle (deg.)    |               24.6               |
            % +-----------------------------+----------------------------------+
            % 
            % Assumptions:
            % - Isotropic bearings;
            %
            % Results:
            %     Reference_5p    Calculated     Rel_Diff      has_Problem
            %     ____________    __________    ___________    ___________
            % 
            %             0              0              NaN       "no"    
            %           710         710.05       -0.0063653       "no"    
            %           710         710.05       -0.0063653       "no"    
            %          1072           1072        -0.001156       "no"    
            %          1072           1072        -0.001156       "no"    
            %        1567.4         1567.4       0.00014933       "no"    
            %        1808.2         1808.2        0.0010573       "no"    
            %        1808.2         1808.2        0.0010573       "no"    
            %        1888.1         1888.1        0.0012111       "no"    
            %        1888.1         1888.1        0.0012111       "no"    
            %        2006.1         2006.1         0.001073       "no"    
            %        2425.3         2425.3       0.00019097       "no"    
            %        2425.3         2425.3       0.00019097       "no"    
            %        2614.8         2614.8       0.00016489       "no"    
            %        5963.8         5963.8       0.00067035       "no"    
            %        5963.8         5963.8       0.00067035       "no"    
            %        6981.7         6981.7      -7.4008e-05       "no"    
            %        6981.7         6981.7      -7.4008e-05       "no"    
            %        7382.4         7382.4       0.00027184       "no"    
            %        7382.4         7382.4       0.00027184       "no"    
            %        8065.4         8065.4       0.00012466       "no"    
            %         11172          11172       -7.746e-05       "no"    
            %         11172          11172       -7.746e-05       "no"    
            %         14253          14253       5.3642e-05       "no"    
            % 
            % Test 01
            % 	Performed on: 24-Nov-2020.
            % 	Number of problematic resonances: 0.
            %

            m_s =  0.4;    m_r =   2.35;   m_c =   5.43;    m_p =   0.66; % [kg]
            d_s = 77.4;    d_r = 275.0;    d_c = 176.8;     d_p = 100.3;  % [mm]
            
            r_s = d_s*1.0e-3/2; % [m]
            r_r = d_r*1.0e-3/2;
            r_c = d_c*1.0e-3/2;
            r_p = d_p*1.0e-3/2;
            
            J_s = 0.39*r_s^2; % [kg-m^2]
            J_r = 3.0 *r_r^2;
            J_c = 6.29*r_c^2;
            J_p = 0.61*r_p^2;
            
            k_sp = 5.0e8;    % k_rp = k_sp; % [N/m]
            k_s  = 1.0e8;   k_r  = k_s;     k_c  = k_s;     k_p  = k_s;  % [N/m]
            k_su = 0.0;     k_ru = 1.0e9;   k_cu = k_su;    k_pu = k_p; % [N/m]
            
            k_st = k_su*r_s^2; % [N-m/rad]
            k_rt = k_ru*r_r^2;
            k_ct = k_cu*r_c^2;
            k_pt = k_pu*r_p^2;
            
            alpha_s = 24.6;     alpha_r = alpha_s;
            
            switch(num_planet)
                case 3
                    fn_ref = sort([0.0 1475.7 1930.3 2658.3 7462.8 11775.3 ...
                                repmat([743.2 1102.4 1896.0 2276.4  6986.3  9647.9], 1, 2)])';
                case 4
                    fn_ref = sort([0.0 1536.6 1970.6 2625.7 7773.6 13071.1 ...
                               repmat([ 727.0 1091.0 1892.8 2342.5  7189.9 10437.6], 1, 2) ...
                                       1808.2 5963.8 6981.7])';
                case 5
                    fn_ref = sort([0.0 1567.4 2006.1 2614.8 8065.4 14253.1 ...
                               repmat([ 710.0 1072.0 1888.1 2425.3  7382.4 11172.3], 1, 2) ...
                               repmat([1808.2 5963.8 6981.7], 1, 2)])';
                otherwise
                    error('Range is between 3 and 5.');
            end
            
            % Defining test cases for:
            % Gear_Set (main):
            import matlab.mock.TestCase;
            gset_TC = TestCase.forInteractiveUse;
            
            % Carrier (aux.):
            carrier_TC = TestCase.forInteractiveUse;
            
            % Creating mock objects and behaviors for:
            % Gear_Set (main):
            [gset_mock, gset_behavior] = createMock(gset_TC, 'AddedProperties',{'configuration', 'N_p', ...
                'alpha_n', 'mass', 'd_b', 'a_w', 'J_x', 'carrier', 'bearing'}, ...
                'AddedMethods', {'sub_set'});
            
            % Carrier (aux.):
            [carrier_mock, carrier_behavior] = createMock(carrier_TC, 'AddedProperties',{'mass', 'J_x'});
            
            % Assigning outputs for:
            % Carrier properties:
            carrier_TC.assignOutputsWhen(get(carrier_behavior.mass), m_c);
            carrier_TC.assignOutputsWhen(get(carrier_behavior.J_x), J_c);
            
            % Gear_Set properties:
            mass = [m_s m_p m_r];
            d    = [d_s d_p d_r];
            J_x  = [J_s J_p J_r];

            brg_s = Bearing('K_y',     k_s, 'K_z',     k_s, 'K_alpha',     k_st, 'name', 'sun'    );
            brg_r = Bearing('K_y',     k_r, 'K_z',     k_r, 'K_alpha',     k_rt, 'name', 'ring'   );
            brg_c = Bearing('K_y', 0.5*k_c, 'K_z', 0.5*k_c, 'K_alpha', 0.5*k_ct, 'name', 'carrier');
            brg_p = Bearing('K_y', 0.5*k_p, 'K_z', 0.5*k_p, 'K_alpha', 0.5*k_pt, 'name', 'planet' );

            brg = [brg_p, brg_p, ...
                   brg_c, brg_c, ...
                   brg_s, ...
                   brg_r];
            
            tmp.k_mesh = k_sp;
            
            gset_TC.assignOutputsWhen(get(gset_behavior.a_w), d_c/2);
            gset_TC.assignOutputsWhen(get(gset_behavior.configuration), 'planetary');
            gset_TC.assignOutputsWhen(get(gset_behavior.N_p), num_planet);
            gset_TC.assignOutputsWhen(get(gset_behavior.alpha_n), alpha_r);
            gset_TC.assignOutputsWhen(get(gset_behavior.mass), mass);
            gset_TC.assignOutputsWhen(get(gset_behavior.d_b), d);
            gset_TC.assignOutputsWhen(get(gset_behavior.J_x), J_x);
            gset_TC.assignOutputsWhen(get(gset_behavior.carrier), carrier_mock);
            gset_TC.assignOutputsWhen(get(gset_behavior.bearing), brg);
            gset_TC.assignOutputsWhen(gset_behavior.sub_set('sun_planet'), tmp);
            gset_TC.assignOutputsWhen(gset_behavior.sub_set('planet_ring'), tmp);
            
            MM = Lin_Parker_99.stage_inertia_matrix(gset_mock);
            [Km, Kb] = Lin_Parker_99.stage_stiffness_matrix(gset_mock);
            KK = Kb + Km;
            
            dyn_for_TC = TestCase.forInteractiveUse;
            id = ["ISO_6336:KV", ...
                  "ISO_6336:SF", ...
                  "ISO_6336:KS", ...
                  "Dynamic_Formulation:imag", ...
                  "Dynamic_Formulation:RB"];
            
            for idx = 1:length(id)
                warning("off", id(idx));
            end
            
            [dyn_for_mock, ~] = createMock(dyn_for_TC, ?Dynamic_Formulation, ...
                                                       'ConstructorInputs', {NREL_5MW('M_R' , 0.0, 'J_R' , 0.0, ...
                                                                                      'M_G' , 0.0, 'J_G' , 0.0)});
            dyn_for_mock.K_m = Km;
            dyn_for_mock.K_b = Kb;
            dyn_for_mock.M = MM;
            dyn_for_mock.A = [zeros(length(MM)),   eye(length(MM));
                                       -MM\KK  , zeros(length(MM))]; 
            
            fn_test = dyn_for_mock.modal_analysis();
            fn_test = sort(fn_test);
            
            for idx = 1:length(id)
                warning("on", id(idx));
            end
            
            ref_str = sprintf('Reference_%dp', num_planet);

			diff_fn = 100*(fn_ref - fn_test)./fn_ref;
            has_problem = repmat("no", length(fn_ref), 1);
            
            idx = abs(diff_fn) > 1.0;
            has_problem(idx) = "[YES]";
            table(fn_ref, fn_test, diff_fn, has_problem, 'variableNames', ...
                    {ref_str, 'Calculated', 'Rel_Diff', 'has_Problem'})
            fprintf('Test 01\n');
            fprintf('\tPerformed on: %s.\n', datetime('today'));
            fprintf('\tNumber of problematic resonances: %d.\n', sum(idx));
            flag = any(idx);
        end
        
        function flag = test_02(num_planet)
            %TEST_02 Testing the code using the parameters from [3], 
            % Table 1, shown below. Reference values for resonances were 
            % taken from [3], Table 2, considering Omega_c = 0.
            %
            % +----------------------------------+-------+-------+---------+--------+
            % |             Property             |  Sun  |  Ring | Carrier | Planet |
            % +----------------------------------+-------+-------+---------+--------+
            % |             Mass (kg)            |  3.00 |  7.64 |  12.00  |  1.86  |
            % +----------------------------------+-------+-------+---------+--------+
            % |            I/r^2 (kg)            |  1.75 |  8.09 |   6.80  |  1.25  |
            % +----------------------------------+-------+-------+---------+--------+
            % |        Base diameter (mm)        | 55.75 | 109.7 |   88.6  |   27   |
            % +----------------------------------+-------+-------+---------+--------+
            % |       Mesh stiffness (N/um)      |                100               |
            % +----------------------------------+-------+-------+---------+--------+
            % | Elastic support stiffness (N/um) |   20  |  100  |    50   |   10   |
            % +----------------------------------+-------+-------+---------+--------+
            % |    Torsional stiffness (N/um)    |   20  |  100  |   100   |    ?   |
            % +----------------------------------+-------+-------+---------+--------+
            % |       Pressure angle (deg.)      |               20.0               |
            % +----------------------------------+----------------------------------+
            %
            % Assumptions:
            % - Isotropic bearing stiffness;
            % - Sun-planet mesh stiffnesses are equal for all planets;
            % - Ring-planet mesh stiffnesses are equal for all planets;
            % - The sun-bearing torsional stiffness is not mentioned and is
            % assumed to be zero.
            %
            % Results:
            %     Reference_5p    Calculated     Rel_Diff     has_Problem
            %     ____________    __________    __________    ___________
            %
            %        258.03         258.04      -0.0062309       "no"
            %        258.03         258.04      -0.0062309       "no"
            %        324.82         324.84      -0.0045065       "no"
            %         354.9         354.91      -0.0046635       "no"
            %         354.9         354.91      -0.0046635       "no"
            %        363.98         363.99      -0.0040313       "no"
            %        363.98         363.99      -0.0040313       "no"
            %        376.04         375.94        0.028082       "no"
            %        472.73         472.53        0.042164       "no"
            %        495.98         496.05       -0.013471       "no"
            %        495.98         496.05       -0.013471       "no"
            %        552.07         552.14       -0.012343       "no"
            %        552.07         552.14       -0.012343       "no"
            %        771.28         771.18        0.013105       "no"
            %        1594.2         1594.1       0.0052395       "no"
            %        1594.2         1594.1       0.0052395       "no"
            %        1962.1         1962.1       0.0043046       "no"
            %        1962.1         1962.1       0.0043046       "no"
            %        2093.1         2093.2       -0.004691       "no"
            %        2093.1         2093.2       -0.004691       "no"
            %        2225.3         2225.4      -0.0053724       "no"
            %        2469.9         2469.9      0.00062585       "no"
            %        2469.9         2469.9      0.00062585       "no"
            %        3316.1           3316       0.0036004       "no"
            %
            % Test 02
            % 	Performed on: 24-Nov-2020.
            % 	Number of problematic resonances: 0.
            %

            switch(num_planet)
                case 3
                    fn_ref = sort([0.9296 1.090 1.358 1.931 5.566 7.741 ...
                           repmat([0.7432 1.018 1.265 1.487 4.976 6.282], 1, 2)]);
                case 4
                    fn_ref = sort([0.9077 1.046 1.317 2.013 5.809 8.383 ...
                           repmat([0.7199 1.001 1.317 1.481 5.153 6.488], 1, 2) ...
                                   0.9617 4.320 5.672]);
                case 5
                    fn_ref = sort([0.8802 1.019  1.281 2.090 6.030 8.986 ...
                           repmat([0.6992 0.9863 1.344 1.496 5.317 6.693], 1, 2) ...
                           repmat([0.9617 4.320  5.672], 1, 2)]);
                otherwise
                    error('Range is between 3 and 5.');
            end
            
            m_s =  3.0;      m_r =   7.64;    m_c =  12.0;    m_p =  1.86; % [kg]
            d_s = 55.75;     d_r = 109.7;     d_c =  88.6;    d_p = 27.0;  % [mm]
            
            r_s = d_s*1.0e-3/2; % [m]
            r_r = d_r*1.0e-3/2;
            r_c = d_c*1.0e-3/2;
            r_p = d_p*1.0e-3/2;
            
            J_s = 1.75*r_s^2; % [kg-m^2]
            J_r = 8.09*r_r^2;
            J_c = 6.80*r_c^2;
            J_p = 1.25*r_p^2;
            
            k_sp = 100.0e6;    % k_rp = k_sp; % [N/m]
            k_s  = 20.0e6;  k_r  = 100.0e6;	k_c  = 50.0e6;	k_p  = 10.0e6;
            k_su = 20.0e6;  k_ru = 100.0e6; k_cu = k_ru;    k_pu = k_p;
            
            k_st = k_su*r_s^2; % [N-m/rad]
            k_rt = k_ru*r_r^2;
            k_ct = k_cu*r_c^2;
            k_pt = k_pu*r_p^2;
            
            omega = sqrt(k_p/m_p)/(2.0*pi);
            fn_ref = fn_ref'*omega;
            
            alpha = 20.0;
            
            % Defining test cases for:
            % Gear_Set (main):
            import matlab.mock.TestCase;
            gset_TC = TestCase.forInteractiveUse;
            
            % Carrier (aux.):
            carrier_TC = TestCase.forInteractiveUse;
            
            % Creating mock objects and behaviors for:
            % Gear_Set (main):
            [gset_mock, gset_behavior] = createMock(gset_TC, 'AddedProperties',{'configuration', 'N_p', ...
                'alpha_n', 'mass', 'd_b', 'a_w', 'J_x', 'carrier', 'bearing', 'output_shaft'}, ...
                'AddedMethods', {'sub_set'});
            
            % Carrier (aux.):
            [carrier_mock, carrier_behavior] = createMock(carrier_TC, 'AddedProperties',{'mass', 'J_x'});
            
            % Assigning outputs for:
            % Carrier properties:
            carrier_TC.assignOutputsWhen(get(carrier_behavior.mass), m_c);
            carrier_TC.assignOutputsWhen(get(carrier_behavior.J_x), J_c);
            
            % Gear_Set properties:
            mass = [m_s m_p m_r];
            d    = [d_s d_p d_r];
            J_x  = [J_s J_p J_r];
            
            brg_s = Bearing('K_y',     k_s, 'K_z',     k_s, 'K_alpha',     k_st, 'name', 'sun'    );
            brg_r = Bearing('K_y',     k_r, 'K_z',     k_r, 'K_alpha',     k_rt, 'name', 'ring'   );
            brg_c = Bearing('K_y', 0.5*k_c, 'K_z', 0.5*k_c, 'K_alpha', 0.5*k_ct, 'name', 'carrier');
            brg_p = Bearing('K_y', 0.5*k_p, 'K_z', 0.5*k_p, 'K_alpha', 0.5*k_pt, 'name', 'planet' );

            brg = [brg_p, brg_p, ...
                   brg_c, brg_c, ...
                   brg_s, ...
                   brg_r];
            
            mesh_stiff.k_mesh = k_sp;
            
            gset_TC.assignOutputsWhen(get(gset_behavior.a_w), d_c/2);
            gset_TC.assignOutputsWhen(get(gset_behavior.configuration), 'planetary');
            gset_TC.assignOutputsWhen(get(gset_behavior.N_p), num_planet);
            gset_TC.assignOutputsWhen(get(gset_behavior.alpha_n), alpha);
            gset_TC.assignOutputsWhen(get(gset_behavior.mass), mass);
            gset_TC.assignOutputsWhen(get(gset_behavior.d_b), d);
            gset_TC.assignOutputsWhen(get(gset_behavior.J_x), J_x);
            gset_TC.assignOutputsWhen(get(gset_behavior.carrier), carrier_mock);
            gset_TC.assignOutputsWhen(get(gset_behavior.bearing), brg);
            gset_TC.assignOutputsWhen(gset_behavior.sub_set('sun_planet'), mesh_stiff);
            gset_TC.assignOutputsWhen(gset_behavior.sub_set('planet_ring'), mesh_stiff);
            
            MM       = Lin_Parker_99.stage_inertia_matrix(gset_mock);
            [Km, Kb] = Lin_Parker_99.stage_stiffness_matrix(gset_mock);
            KK = Kb + Km;
            
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
            dyn_for_mock.K_b = Kb;
            dyn_for_mock.K_m = Km;
            dyn_for_mock.M = MM;
            dyn_for_mock.A = [zeros(length(MM)),   eye(length(MM));
                                       -MM\KK  , zeros(length(MM))]; 
            
            fn_test = dyn_for_mock.modal_analysis();
            fn_test = sort(fn_test);
            
            for idx = 1:length(id)
                warning("on", id(idx));
            end
            
            ref_str = sprintf('Reference_%dp', num_planet);

			diff_fn = 100*(fn_ref - fn_test)./fn_ref;
            has_problem = repmat("no", length(fn_ref), 1);
            
            idx = abs(diff_fn) > 1.0;
            has_problem(idx) = "[YES]";
            table(fn_ref, fn_test, diff_fn, has_problem, 'variableNames', ...
                    {ref_str, 'Calculated', 'Rel_Diff', 'has_Problem'})
            fprintf('Test 02\n');
            fprintf('\tPerformed on: %s.\n', datetime('today'));
            fprintf('\tNumber of problematic resonances: %d.\n', sum(idx));
            flag = any(idx);
        end
        
    end
    
    %% get methods:
    methods
        function val = get.KK(obj)
            val = @(Om)(obj.K_Bearing + obj.K_mesh - obj.K_Omega*Om.^2);
        end
    end
end
