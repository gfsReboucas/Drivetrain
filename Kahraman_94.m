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
            
            obj.n_DOF = obj.calculate_DOF();
            obj.M = obj.inertia_matrix();
            obj.K = obj.stiffness_matrix();
            
            obj.load      = zeros(obj.n_DOF(end), 1);
            obj.load(1)   = 1.0;
            obj.load(end) = -obj.load(1);
        end
        
        %% Calculation:
        function KK = stiffness_matrix(obj)
            N = obj.n_DOF(end);
            KK = zeros(N, N);
            
            K_tmp = obj.main_shaft.stiffness_matrix("torsional");
            KK(1:2, 1:2) = KK(1:2, 1:2) + K_tmp;
            
            for idx = 1:obj.N_stage
                K_tmp = Kahraman_94.stage_stiffness_matrix(obj.stage(idx));
                jdx = obj.n_DOF(idx);
                kdx = jdx:(jdx + length(K_tmp) - 1);
                
                KK(kdx, kdx) = KK(kdx, kdx) + K_tmp;
            end
        end
        
        function MM = inertia_matrix(obj)
            N = obj.n_DOF(end);
            MM = zeros(N, N);
            
            MM(1, 1) = obj.J_Rotor;     MM(end, end) = obj.J_Gen;
            
            M_tmp = obj.main_shaft.inertia_matrix("torsional");
            MM(1:2, 1:2) = MM(1:2, 1:2) + M_tmp;
            
            for idx = 1:obj.N_stage
                M_tmp = Kahraman_94.stage_inertia_matrix(obj.stage(idx));
                jdx = obj.n_DOF(idx);
                kdx = jdx:(jdx + length(M_tmp) - 1);
                
                MM(kdx, kdx) = MM(kdx, kdx) + M_tmp;
            end
        end
        
        function nn = calculate_DOF(obj)
            nn = ones(obj.N_stage + 1, 1)*2.0;
            
            for idx = 1:obj.N_stage
                if(strcmp(obj.stage(idx).configuration, "parallel"))
                    tmp = obj.stage(idx).N_p + 1;
                elseif(strcmp(obj.stage(idx).configuration, "planetary"))
                    tmp = obj.stage(idx).N_p + 2;
                end
                
                nn(idx + 1) = nn(idx) + tmp;
            end
        end
    end
    
    methods(Static)
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

                k_sp = abs(sun_pla.k_mesh);
                k_rp = abs(pla_rng.k_mesh); % got a negative k_mesh for DTU_10MW

                aw  =  stage_idx.a_w *1.0e-3;
                r_s = (stage_idx.d(1)*1.0e-3)/2.0;
                r_p = (stage_idx.d(2)*1.0e-3)/2.0;
                
                KK(    1,     1) =  stage_idx.N_p*(k_rp + k_sp)*aw^2;
                KK(    1, n - 1) = -stage_idx.N_p*  r_s * k_sp *aw;
                KK(n - 1,     1) = KK(1, n - 1);
                KK(n - 1, n - 1) =  stage_idx.N_p*k_sp*r_s^2;
                
                for jdx = 2:(n - 2)
                    KK(    1,   jdx) = aw*r_p*(k_rp - k_sp);
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
                
                aw  =  stage_idx.a_w *1.0e-3;
                r_s = (stage_idx.d(1)*1.0e-3)/2.0;
                r_p = (stage_idx.d(2)*1.0e-3)/2.0;
                
                range = 1:(n - 1);
                MM(range, range) = diag([m_c* aw^2, ...
                                         m_p*r_p^2*ones(1, stage_idx.N_p), ...
                                         m_s*r_s^2]);
            end
            
            range = (n - 1):n;
            
            MM(range, range) = MM(range, range) + ...
                stage_idx.output_shaft.inertia_matrix('torsional');
        end
        
    end
end

