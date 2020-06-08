classdef Lin_Parker_99 < Dynamic_Formulation
    %KAHRAMAN_94 Calculates the inertia and stiffness matrix of a
    % multi-stage Drivetrain object according to:
    % [1] J. Lin and R. Parker, 'Analytical Characterization of the
    % Unique Properties of Planetary Gear Free Vibration', Journal
    % of Vibration and Acoustics, vol. 121, no. 3, pp. 316-321,
    % 1999. https://doi.org/10.1115/1.2893982
    %
    properties
        K_bearing; % Bearing stiffness
        K_mesh;    % Mesh stiffness
        K_Omega;   % Centripetal stiffness
        G;         % Gyroscopic matrix
    end
    
    properties(Dependent)
        KK;
    end
    
    methods
        function obj = Lin_Parker_99(DT)
            obj@Dynamic_Formulation(DT);
            
            obj.n_DOF = obj.calculate_DOF();
            [obj.M, obj.G] = obj.inertia_matrix();
            [obj.K_bearing, obj.K_mesh, obj.K_Omega] = obj.stiffness_matrix();
            obj.K = obj.K_bearing + obj.K_mesh;
            
            obj.load      = zeros(obj.n_DOF(end), 1);
            obj.load(1:3) = 1.0;
            obj.load(end-2:end) = -obj.load(1:3);
        end
        
        %% Calculation:
        function [Kb, Km, KO] = stiffness_matrix(obj)
            N = obj.n_DOF(end);
            Kb = zeros(N, N);
            Km = zeros(N, N);
            KO = zeros(N, N);
            
            K_tmp = obj.main_shaft.stiffness_matrix('Lin_Parker_99');
            
            range = 1:6;
            Kb(range, range) = Kb(range, range) + K_tmp;

            N = obj.n_DOF;
            
            for idx = 1:obj.N_stage
                [Kbt, Kmt, KOt] = obj.stage_stiffness_matrix(idx);
                
                range = (N(idx) - 2):(N(idx + 1));

                Kb(range, range) = Kb(range, range) + Kbt;
                Km(range, range) = Km(range, range) + Kmt;
                KO(range, range) = KO(range, range) + KOt;
            end
        end
        
        function [Kb, Km, KO] = stage_stiffness_matrix(obj, idx)
            stage_idx = obj.stage(idx);
            alpha_n = stage_idx.alpha_n;
            
            N1    = stage_idx.N_p - 1;
            psi   = (0:N1)*(360.0/stage_idx.N_p);
            psi_s = psi - alpha_n;
            
            % Bearing stiffness sub-matrix:
            K_b_ = @(x, y)diag([x, y, 0]);
            
            % sun-sun mesh-stiffness matrix:
            K_s1 = @(i)[                sind(psi_s(i))^2, -cosd(psi_s(i))*sind(psi_s(i)), -sind(psi_s(i));
                        -cosd(psi_s(i))*sind(psi_s(i))  ,  cosd(psi_s(i))^2             ,  cosd(psi_s(i));
                        -               sind(psi_s(i))  ,  cosd(psi_s(i))               ,         1     ];
            
            % sun-planet mesh-stiffness matrix:
            K_s2 = @(i)[ sind(psi_s(i))*sind(alpha_n),  sind(psi_s(i))*cosd(alpha_n), -sind(psi_s(i));
                        -cosd(psi_s(i))*sind(alpha_n), -cosd(psi_s(i))*cosd(alpha_n),  cosd(psi_s(i));
                        -               sind(alpha_n), -               cosd(alpha_n),         1     ];
                                 
            % planet-planet (?) mesh-stiffness matrix:
            K_s3 = [ sind(alpha_n)^2            ,  sind(alpha_n)*cosd(alpha_n)  , -sind(alpha_n);
                     sind(alpha_n)*cosd(alpha_n),                cosd(alpha_n)^2, -cosd(alpha_n);
                    -sind(alpha_n)              , -              cosd(alpha_n)  ,        1     ];
            
            % carrier-carrier bearing stiffness matrix:
            K_c1 = @(i)[ 1           , 0           , -sind(psi(i));
                         0           , 1           ,  cosd(psi(i));
                        -sind(psi(i)), cosd(psi(i)),        1    ];

            % carrier-planet bearing stiffness matrix:
            K_c2 = @(i)[-cosd(psi(i)),  sind(psi(i)), 0;
                        -sind(psi(i)), -cosd(psi(i)), 0;
                         0           , -1           , 0];
            
            % (?)
            K_c3 = @(x, y)K_b_(x, y);
            
            % (?)
            K_r3 = [ sind(alpha_n)^2            , -sind(alpha_n)*cosd(alpha_n)  , -sind(alpha_n);
                    -sind(alpha_n)*cosd(alpha_n),                cosd(alpha_n)^2,  cosd(alpha_n);
                    -sind(alpha_n)              ,                cosd(alpha_n)  ,        1     ];
            
            if(strcmp(stage_idx.configuration, 'parallel'))
                n = 9;
                Kb = zeros(n, n);
                Km = zeros(n, n);
                KO = zeros(n, n);
                R  = zeros(n, n);
                
                range = 1:(n - 3);
                
                % Bearing component:
                b_p  = stage_idx.bearing(4:end);
                b_p  = b_p.parallel_association();
                k_px = b_p.K_y;
                k_py = b_p.K_z;
                
                Kb(range, range) = diag([zeros(1,3) k_px k_py 0]);
                
                % Mesh component:
                b_w  = stage_idx.bearing(1:3);
                b_w  = b_w.parallel_association();
                k_wx = b_w.K_y;
                k_wy = b_w.K_z;
                
                k = stage_idx.k_mesh;
                
                Km(range, range) = [k*K_s3 + K_c3(k_wx, k_wy), k*K_s2(1) ;
                                    k*K_s2(1)                , k*K_s1(1)];
                
                % Centripetal component:
                m_1 = stage_idx.mass(1);
                m_2 = stage_idx.mass(2);
                KO(range, range) = diag([m_2, m_2, 0, ... % wheel
                                         m_1, m_1, 0]);   % pinion

                % Torsional to translational:
                r_p = (stage_idx.d(1)*1.0e-3)/2.0;
                r_w = (stage_idx.d(2)*1.0e-3)/2.0;
                
                R(range, range) = diag([1.0, 1.0, r_w, ... % wheel
                                        1.0, 1.0, r_p]);   % pinion
            elseif(strcmp(stage_idx.configuration, 'planetary'))
                n = 3*(stage_idx.N_p + 3);
                Kb = zeros(n, n);
                Km = zeros(n, n);
                KO = zeros(n, n);
                R  = zeros(n, n);
                
                np = 3*stage_idx.N_p;
                range = 1:(n - 3);
                
                % Bearing component:
                b_c = stage_idx.bearing(3:end);
                b_c = b_c.parallel_association();
                k_cx = b_c.K_y;
                k_cy = b_c.K_z;
                
                Kb(range, range) = diag([k_cx, k_cy, 0.0, ... % carrier
                                         zeros(1, np),    ... % planet
                                         zeros(1, 3)]);       % sun
                
                % Mesh component:
                b_p = stage_idx.bearing(1:2);
                b_p = b_p.parallel_association();
                
                k_px = b_p.K_y;
                k_py = b_p.K_z;
                
                sun_pla = Gear_Set('configuration', 'parallel'               , ...
                                   'm_n'          , stage_idx.m_n            , ...
                                   'alpha_n'      , stage_idx.alpha_n        , ...
                                   'z'            , stage_idx.z(1:2)         , ...
                                   'b'            , stage_idx.b              , ...
                                   'beta'         , stage_idx.beta           , ...
                                   'x'            , stage_idx.x(1:2)         , ...
                                   'k'            , stage_idx.k(1:2)         , ...
                                   'bore_ratio'   , stage_idx.bore_ratio(1:2), ...
                                   'N_p'          , 1                        , ...
                                   'a_w'          , stage_idx.a_w            , ...
                                   'rack_type'    , stage_idx.type           , ...
                                   'bearing'      , stage_idx.bearing        , ...
                                   'shaft'        , stage_idx.output_shaft);

                pla_rng = Gear_Set('configuration', 'parallel'               , ...
                                   'm_n'          , stage_idx.m_n            , ...
                                   'alpha_n'      , stage_idx.alpha_n        , ...
                                   'z'            , stage_idx.z(2:3)         , ...
                                   'b'            , stage_idx.b              , ...
                                   'beta'         , stage_idx.beta           , ...
                                   'x'            , stage_idx.x(2:3)         , ...
                                   'k'            , stage_idx.k(2:3)         , ...
                                   'bore_ratio'   , stage_idx.bore_ratio(2:3), ...
                                   'N_p'          , 1                        , ...
                                   'a_w'          , stage_idx.a_w            , ...
                                   'rack_type'    , stage_idx.type           , ...
                                   'bearing'      , stage_idx.bearing        , ...
                                   'shaft'        , stage_idx.output_shaft);

                % OBS.: got a negative k_mesh for DTU_10MW planet-ring
                k_sp = sun_pla.k_mesh;
                k_pr = pla_rng.k_mesh;
                
                K_c = zeros(3, np);
                K_s = K_c';
                sum_Kc = 0;
                sum_Ks = 0;
                for idx = 1:stage_idx.N_p
                    rng = (-2:0) + 3*idx;
                    K_c(:, rng) = K_c2(idx);
                    K_s(rng, :) = K_s2(idx);
                    
                    sum_Kc = sum_Kc + K_c1(idx);
                    sum_Ks = sum_Ks + K_s1(idx);
                end
                
                K_c    = K_b_(k_px, k_py)*K_c;
                K_s    =            k_sp *K_s;
                sum_Kc = K_b_(k_px, k_py)*sum_Kc;
                sum_Ks =            k_sp *sum_Ks;

                K_pp = K_c3(k_cx, k_cy) + k_pr*K_r3 + k_sp*K_s3;
                KppCell = repmat({K_pp}, 1, stage_idx.N_p);
                Kpp = blkdiag(KppCell{:});
                
                np3 = np + 3;
                Km(1:3  , 4:np3) = K_c;
                Km(4:np3, (1:3) + np+3) = K_s;
                Km(range, range) = Km(range, range) + Km(range, range)';
                Km(range, range) = Km(range, range) + blkdiag(sum_Kc, ... % carrier
                                                              Kpp, ...    % planet
                                                              sum_Ks);    % sun

                % Centripetal component:
                m_c = stage_idx.carrier.mass;
                m_s = stage_idx.mass(1);
                m_p = stage_idx.mass(2);
                                                          
                KO(range, range) = diag([m_c m_c 0.0, ...
                                        [m_p m_p 0.0]*repmat(eye(3), 1, stage_idx.N_p), ...
                                         m_s m_s 0.0]);
                
                % Torsional to translational:
                r_c =  stage_idx.a_w *1.0e-3;
                r_s = (stage_idx.d(1)*1.0e-3)/2.0;
                r_p = (stage_idx.d(2)*1.0e-3)/2.0;
                
                R(range, range) = diag([1.0, 1.0, r_c, ...
                                       [1.0, 1.0, r_p]*repmat(eye(3), 1, stage_idx.N_p), ...
                                        1.0, 1.0, r_s]);
            end
            
            range = (n - 2):n;
            R(range, range) = eye(3);
            
            range = (n - 5):n;
            Kb(range, range) = Kb(range, range) + ...
                stage_idx.output_shaft.stiffness_matrix('Lin_Parker_99');
            
            % Torsional to translational:
            Kb = R' * Kb * R;
            Km = R' * Km * R;
            KO = R' * KO * R;
            
            % making it symmetric, disregarding spurious values:
            Kb = (Kb + Kb')./2.0;
            Km = (Km + Km')./2.0;
            KO = (KO + KO')./2.0;
        end
        
        function [MM, GG] = inertia_matrix(obj)
            N = obj.n_DOF(end);
            MM = zeros(N, N);
            GG = zeros(N, N);
            
            m_r = obj.m_Rotor;            m_g = obj.m_Gen;
            J_r = obj.J_Rotor;            J_g = obj.J_Gen;
            
            range = 1:3;
            MM(range, range) = diag([m_r, m_r, J_r]);
            
            range = (N-2):N;
            MM(range, range) = diag([m_g, m_g, J_g]);
            
            M_tmp = obj.main_shaft.inertia_matrix('Lin_Parker_99');
            
            range = 1:6;
            MM(range, range) = MM(range, range) + M_tmp;
            
            N = obj.n_DOF;
            
            for idx = 1:obj.N_stage
                [M_tmp, G_tmp] = obj.stage_inertia_matrix(idx);
                
                range = (N(idx) - 2):(N(idx + 1));
                
                MM(range, range) = MM(range, range) + M_tmp;
                GG(range, range) = GG(range, range) + G_tmp;
            end
        end
        
        function [MM, GG] = stage_inertia_matrix(obj, idx)
            stage_idx = obj.stage(idx);
            G_ = @(x)[0.0  , -2.0*x, 0.0;
                      2.0*x,  0.0  , 0.0;
                      0.0  ,  0.0  , 0.0];
            
            if(strcmp(stage_idx.configuration, 'parallel'))
                n = 9;
                
                m_p = stage_idx.mass(1);
                m_w = stage_idx.mass(2);
                
                r_p = (stage_idx.d(1)*1.0e-3)/2.0;
                r_w = (stage_idx.d(2)*1.0e-3)/2.0;
                
                J_p = stage_idx.J_x(1)/r_p^2;
                J_w = stage_idx.J_x(2)/r_w^2;
                
                d = [m_w, m_w, J_w, ... % wheel
                     m_p, m_p, J_p, ... % pinion
                     0.0, 0.0, 0.0];    % shaft
                
                % Gyroscopic matrix:
                GG = blkdiag(G_(m_w), G_(m_p), zeros(3));
                
                % Torsional to translational:
                R = diag([1.0, 1.0, r_w, ... % wheel
                          1.0, 1.0, r_p, ... % pinion
                          1.0, 1.0, 1.0]);   % shaft

            elseif(strcmp(stage_idx.configuration, 'planetary'))
                n = 3*(stage_idx.N_p + 3);
                
                m_s = stage_idx.mass(1);
                m_p = stage_idx.mass(2);
                m_c = stage_idx.carrier.mass;
                
                r_s = (stage_idx.d(1)*1.0e-3)/2.0;
                r_p = (stage_idx.d(2)*1.0e-3)/2.0;
                r_c =  stage_idx.a_w *1.0e-3;
                
                J_s = stage_idx.J_x(1)/r_s^2;
                J_p = stage_idx.J_x(2)/r_p^2;
                J_c = stage_idx.carrier.J_x/r_c^2;
                
                d = [m_c, m_c, J_c, ...
                    [m_p, m_p, J_p]*repmat(eye(3), 1, stage_idx.N_p), ...
                     m_s, m_s, J_s, ...
                     0.0, 0.0, 0.0];
                
                % Gyroscopic matrix:
                Gp = repmat({G_(m_p)}, 1, stage_idx.N_p);
                Gp = blkdiag(Gp{:});
                
                GG = blkdiag(G_(m_c), ... % carrier
                             Gp     , ... % planet
                             G_(m_s), ... % sun
                             G_(0.0));    % shaft
                
                % Torsional to translational:
                R = diag([1.0, 1.0, r_c, ...
                         [1.0, 1.0, r_p]*repmat(eye(3), 1, stage_idx.N_p), ...
                          1.0, 1.0, r_s, ...
                          1.0, 1.0, 1.0]);
            end
            
            MM = diag(d);
            
            range = (n - 5):n;
            MM(range, range) = MM(range, range) + ...
                stage_idx.output_shaft.inertia_matrix('Lin_Parker_99');
            
            MM = R' * MM * R;
            GG = R' * GG * R;
            
        end
        
        function nn = calculate_DOF(obj)
            nn = ones(obj.N_stage + 1, 1)*6.0;
            
            for idx = 1:obj.N_stage
                if(strcmp(obj.stage(idx).configuration, 'parallel'))
                    tmp = obj.stage(idx).N_p + 1;
                elseif(strcmp(obj.stage(idx).configuration, 'planetary'))
                    tmp = obj.stage(idx).N_p + 2;
                end
                
                nn(idx + 1) = nn(idx) + tmp*3;
            end
        end
    end
    
    %% get methods:
    methods
        function val = get.KK(obj)
            val = @(Om)(obj.K_Bearing + obj.K_mesh - obj.K_Omega*Om.^2);
        end
    end
end
