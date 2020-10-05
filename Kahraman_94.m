classdef Kahraman_94 < Dynamic_Formulation
    %KAHRAMAN_94 Calculates the inertia and stiffness matrices of a
    % multi-stage Drivetrain object according to:
    % [1] A. Kahraman, "Natural Modes of Planetary Gear Trains",
    % Journal of Sound and Vibration, vol. 173, no. 1, pp. 125-130,
    % 1994. https://doi.org/10.1006/jsvi.1994.1222.
    %
    
    methods
        function obj = Kahraman_94(DT)
            obj@Dynamic_Formulation(DT);
        end
        
        %% Calculation:
        function MM = inertia_matrix(obj)
            N = obj.n_DOF(end);
            MM = zeros(N, N);
            
            DT = obj.drive_train;

            MM(1, 1) = DT.J_Rotor;     MM(end, end) = DT.J_Gen;
            
            M_tmp = DT.main_shaft.inertia_matrix("torsional");
            MM(1:2, 1:2) = MM(1:2, 1:2) + M_tmp;
            
            for idx = 1:DT.N_stage
                M_tmp = Kahraman_94.stage_inertia_matrix(DT.stage(idx));
                jdx = obj.n_DOF(idx);
                kdx = jdx:(jdx + length(M_tmp) - 1);
                
                MM(kdx, kdx) = MM(kdx, kdx) + M_tmp;
            end
        end
        
        function KK = stiffness_matrix(obj)
            N = obj.n_DOF(end);
            KK = zeros(N, N);
            
            DT = obj.drive_train;

            K_tmp = DT.main_shaft.stiffness_matrix("torsional");
            KK(1:2, 1:2) = KK(1:2, 1:2) + K_tmp;
            
            for idx = 1:DT.N_stage
                K_tmp = Kahraman_94.stage_stiffness_matrix(DT.stage(idx));
                jdx = obj.n_DOF(idx);
                kdx = jdx:(jdx + length(K_tmp) - 1);
                
                KK(kdx, kdx) = KK(kdx, kdx) + K_tmp;
            end
        end
        
        function DD = damping_matrix(obj)
            N = obj.n_DOF(end);
            DD = zeros(N, N);
            
            DT = obj.drive_train;
            
            D_tmp = DT.main_shaft.damping_matrix("torsional");
            DD(1:2, 1:2) = DD(1:2, 1:2) + D_tmp;
            
            for idx = 1:DT.N_stage
                D_tmp = Kahraman_94.stage_damping_matrix(DT.stage(idx));
                jdx = obj.n_DOF(idx);
                kdx = jdx:(jdx + length(D_tmp) - 1);
                
                DD(kdx, kdx) = DD(kdx, kdx) + D_tmp;
            end
        end
        
        function bb = load_vector(obj)
            bb         = zeros(obj.n_DOF(end), 2);
            bb(1  , 1) = 1.0;
            bb(end, 2) = 1.0;
        end
        
        function nn = calculate_DOF(obj)
            DT = obj.drive_train;
            nn = ones(DT.N_stage + 1, 1)*2.0;
            
            for idx = 1:DT.N_stage
                if(strcmp(DT.stage(idx).configuration, "parallel"))
                    tmp = DT.stage(idx).N_p + 1;
                elseif(strcmp(DT.stage(idx).configuration, "planetary"))
                    tmp = DT.stage(idx).N_p + 2;
                end
                
                nn(idx + 1) = nn(idx) + tmp;
            end
        end
        
    end
    
    methods(Static)
        function MM = stage_inertia_matrix(stage_idx)
            
            if(strcmp(stage_idx.configuration, 'parallel'))
                n = 3;
                MM = zeros(n, n);
                
                m_p = stage_idx.mass(1);
                m_w = stage_idx.mass(2);
                
                r_p = (stage_idx.d(1)*1.0e-3)/2.0;
                r_w = (stage_idx.d(2)*1.0e-3)/2.0;
                
                range = 1:2;
                MM(range, range) = diag([m_w*r_w^2 m_p*r_p^2]);
            elseif(strcmp(stage_idx.configuration, 'planetary'))
                n = stage_idx.N_p + 3;
                MM = zeros(n, n);
                
                m_s = stage_idx.mass(1);
                m_p = stage_idx.mass(2);
                m_c = stage_idx.carrier.mass;
                
                r_c =  stage_idx.a_w *1.0e-3;
                r_s = (stage_idx.d(1)*1.0e-3)/2.0;
                r_p = (stage_idx.d(2)*1.0e-3)/2.0;
                
                range = 1:(n - 1);
                MM(range, range) = diag([m_c*r_c^2, ...
                                         m_p*r_p^2*ones(1, stage_idx.N_p), ...
                                         m_s*r_s^2]);
            end
            
            range = (n - 1):n;
            
            MM(range, range) = MM(range, range) + ...
                stage_idx.output_shaft.inertia_matrix('torsional');
        end
        
        function KK = stage_stiffness_matrix(stage_idx)

            if(strcmp(stage_idx.configuration, 'parallel'))
                n = 3;
                KK = zeros(n, n);
                
                r_p = (stage_idx.d(1)*1.0e-3)/2.0;
                r_w = (stage_idx.d(2)*1.0e-3)/2.0;
                
                k_pw = stage_idx.k_mesh;
                
                range = 1:2;
                KK(range, range) = [r_w ^ 2 * k_pw         r_p     * k_pw * r_w;
                                r_w     * k_pw * r_p   r_p ^ 2 * k_pw];
            elseif(strcmp(stage_idx.configuration, 'planetary'))
                n = stage_idx.N_p + 3;
                KK = zeros(n, n);
                
                % Mesh stiffness:
                % Mesh component:
                sun_pla = stage_idx.sub_set('sun_planet');
                pla_rng = stage_idx.sub_set('planet_ring');

                
                % OBS.: got a negative k_mesh for DTU_10MW planet-ring
                k_sp = sun_pla.k_mesh;
                k_rp = pla_rng.k_mesh;

                r_c =  stage_idx.a_w *1.0e-3;
                r_s = (stage_idx.d(1)*1.0e-3)/2.0;
                r_p = (stage_idx.d(2)*1.0e-3)/2.0;
                
                KK(    1,     1) =  stage_idx.N_p*(k_rp + k_sp)*r_c^2;
                KK(    1, n - 1) = -stage_idx.N_p*  r_s * k_sp *r_c;
                KK(n - 1,     1) = KK(1, n - 1);
                KK(n - 1, n - 1) =  stage_idx.N_p*k_sp*r_s^2;
                
                for jdx = 2:(n - 2)
                    KK(    1,   jdx) = r_c*r_p*(k_rp - k_sp);
                    KK(jdx  ,     1) = KK(1, jdx);
                    KK(jdx  ,   jdx) = (k_rp + k_sp)*r_p^2;
                    KK(n - 1,   jdx) = r_s*r_p*k_sp;
                    KK(jdx  , n - 1) = KK(n - 1, jdx);
                end
            end
            range = (n - 1):n;
            
            KK(range, range) = KK(range, range) + ...
                stage_idx.output_shaft.stiffness_matrix('torsional');
            
        end
        
        function DD = stage_damping_matrix(stage_idx)

            if(strcmp(stage_idx.configuration, 'parallel'))
                n = 3;
                DD = zeros(n, n);
                
                r_p = (stage_idx.d(1)*1.0e-3)/2.0;
                r_w = (stage_idx.d(2)*1.0e-3)/2.0;
                
                d_pw = 500.0e6;
                
                range = 1:2;
                DD(range, range) = [r_w ^ 2 * d_pw         r_p     * d_pw * r_w;
                                    r_w     * d_pw * r_p   r_p ^ 2 * d_pw];
            elseif(strcmp(stage_idx.configuration, 'planetary'))
                n = stage_idx.N_p + 3;
                DD = zeros(n, n);
                
                d_sp = 500.0e6;
                d_rp = 500.0e6;
                
                r_c =  stage_idx.a_w *1.0e-3;
                r_s = (stage_idx.d(1)*1.0e-3)/2.0;
                r_p = (stage_idx.d(2)*1.0e-3)/2.0;
                
                DD(    1,     1) =  stage_idx.N_p*(d_rp + d_sp)*r_c^2;
                DD(    1, n - 1) = -stage_idx.N_p* r_s  * d_sp *r_c;
                DD(n - 1,     1) = DD(1, n - 1);
                DD(n - 1, n - 1) =  stage_idx.N_p*d_sp*r_s^2;
                
                for jdx = 2:(n - 2)
                    DD(    1,   jdx) = r_c*r_p*(d_rp - d_sp);
                    DD(jdx  ,     1) = DD(1, jdx);
                    DD(jdx  ,   jdx) = (d_rp + d_sp)*r_p^2;
                    DD(n - 1,   jdx) = r_s*r_p*d_sp;
                    DD(jdx  , n - 1) = DD(n - 1, jdx);
                end
            end
            range = (n - 1):n;
            
            DD(range, range) = DD(range, range) + ...
                stage_idx.output_shaft.damping_matrix('torsional');
            
        end
        
    end
end

