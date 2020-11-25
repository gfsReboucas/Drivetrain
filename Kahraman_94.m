classdef Kahraman_94 < Dynamic_Formulation
    %KAHRAMAN_94 Calculates the inertia and stiffness matrices of a
    % multi-stage Drivetrain object according to:
    % [1] A. Kahraman, "Natural Modes of Planetary Gear Trains",
    % Journal of Sound and Vibration, vol. 173, no. 1, pp. 125-130,
    % 1994. https://doi.org/10.1006/jsvi.1994.1222.
    %
    
    methods
        function obj = Kahraman_94(DT, varargin)
            default = {'fault_type',  '', ...
                       'fault_val' , 0.0};
            
            default = scaling_factor.process_varargin(default, varargin);
            
            obj@Dynamic_Formulation(DT, 'fault_type', default.fault_type, ...
                                        'fault_val' , default.fault_val);
        end
        
        function desc = explain_DOF(obj)
            N = obj.n_DOF(end);
            desc = cell(2*N, 2);
            jdx = 1;
            desc{jdx, 1} = 'Rotor angular displacement, [rad]';
            desc{jdx, 2} = 'theta_R';
            desc{N  , 1} = 'Generator angular displacement, [rad]';
            desc{N  , 2} = 'theta_G';
            
            DT = obj.drive_train;
            
            for idx = 1:DT.N_stage
                stage_idx = DT.stage(idx);
                if(strcmp(stage_idx.configuration, 'parallel'))
                    jdx = jdx + 1;
                    desc{jdx, 1} = sprintf('Stage %d: Wheel angular displacement, [rad]', idx);
                    desc{jdx, 2} = sprintf('theta_W%d', idx);
                    jdx = jdx + 1;
                    desc{jdx, 1} = sprintf('Stage %d: Pinion angular displacement, [rad]', idx);
                    desc{jdx, 2} = sprintf('theta_P%d', idx);
                elseif(strcmp(stage_idx.configuration, 'planetary'))
                    jdx = jdx + 1;
                    desc{jdx, 1} = sprintf('Stage %d: Carrier angular displacement, [rad]', idx);
                    desc{jdx, 2} = sprintf('theta_c%d', idx);
                    for kdx = 1:stage_idx.N_p
                        desc{jdx + kdx, 1} = sprintf('Stage %d: Planet %d angular displacement, [rad]', idx, kdx);
                        desc{jdx + kdx, 2} = sprintf('theta_p%d%d', idx, kdx);
                    end
                    jdx = jdx + kdx + 1;
                    desc{jdx, 1} = sprintf('Stage %d: Sun angular displacement, [rad]', idx);
                    desc{jdx, 2} = sprintf('theta_s%d', idx);
                end
            end
            
            for idx = 1:N
                desc{N + idx, 1} = strrep(desc{idx, 1}, 'displacement, [rad]', ...
                                                        'speed, [1/min]');
                desc{N + idx, 2} = strrep(desc{idx, 2}, 'theta', ...
                                                        'omega');
            end
        end
        
        %% Calculation:
        function [MM, GG] = inertia_matrix(obj)
            N = obj.n_DOF(end);
            MM = zeros(N);      GG = zeros(N);
            
            DT = obj.drive_train;

            MM(1, 1) = DT.J_Rotor;     MM(end, end) = DT.J_Gen;
            
            M_tmp = DT.main_shaft.inertia_matrix("torsional");
            MM(1:2, 1:2) = MM(1:2, 1:2) + M_tmp;
            
            for idx = 1:DT.N_stage
                if(strcmpi(obj.fault_type, 'M') && (obj.fault_stage == idx))
                    M_tmp = Kahraman_94.stage_faulty_inertia_matrix(DT.stage(idx), ...
                                                                    obj.fault_val(idx));
                else
                    M_tmp = Kahraman_94.stage_inertia_matrix(DT.stage(idx));
                end
                jdx = obj.n_DOF(idx);
                kdx = jdx:(jdx + length(M_tmp) - 1);
                
                MM(kdx, kdx) = MM(kdx, kdx) + M_tmp;
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
        
        function [Km, Kb, KOm] = stiffness_matrix(obj)
            N = obj.n_DOF(end);
            Kb  = zeros(N);            KOm = zeros(N);
            
            DT = obj.drive_train;
            
            % Mesh stiffness matrix:
            Km = zeros(N);
            
            K_tmp = DT.main_shaft.stiffness_matrix("torsional");
            Km(1:2, 1:2) = Km(1:2, 1:2) + K_tmp;
            
            for idx = 1:DT.N_stage
                if(strcmpi(obj.fault_type, 'K') && (obj.fault_stage == idx))
                    K_tmp = Kahraman_94.stage_faulty_stiffness_matrix(DT.stage(idx), ...
                                                                      obj.fault_val(idx));
                else
                    K_tmp = Kahraman_94.stage_stiffness_matrix(DT.stage(idx));
                end
                jdx = obj.n_DOF(idx);
                kdx = jdx:(jdx + length(K_tmp) - 1);
                
                Km(kdx, kdx) = Km(kdx, kdx) + K_tmp;
            end
        end
        
        function cc = centripetal_force_vector(obj)
            cc = zeros(obj.n_DOF(end), 1);
        end
        
        function bb = external_load_vector(obj)
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
        
        function MM = stage_faulty_inertia_matrix(stage_idx, val)
            
            if(strcmp(stage_idx.configuration, 'parallel'))
                n = 3;
                MM = zeros(n, n);
                
                m_p = stage_idx.mass(1);
                m_w = stage_idx.mass(2);
                m_w = m_w*(1.0 - val);
                
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
                
                vec = ones(1, stage_idx.N_p);
                vec(2) = 1.0 - val;
                
                range = 1:(n - 1);
                MM(range, range) = diag([m_c*r_c^2, ...
                                         m_p*r_p^2*vec, ...
                                         m_s*r_s^2]);
            end
            
            range = (n - 1):n;
            
            MM(range, range) = MM(range, range) + ...
                stage_idx.output_shaft.inertia_matrix('torsional');
        end
        
        function DD = stage_damping_matrix(stage_idx)
            
            mesh_damp = 500.0e6;
%             mesh_damp = 0.0;
            if(strcmp(stage_idx.configuration, 'parallel'))
                n = 3;
                DD = zeros(n, n);
                
                r_p = (stage_idx.d(1)*1.0e-3)/2.0;
                r_w = (stage_idx.d(2)*1.0e-3)/2.0;
                
                d_pw = mesh_damp; %*stage_idx.k_mesh;
                
                range = 1:2;
                DD(range, range) = [r_w ^ 2 * d_pw         r_p     * d_pw * r_w;
                                    r_w     * d_pw * r_p   r_p ^ 2 * d_pw];
            elseif(strcmp(stage_idx.configuration, 'planetary'))
                n = stage_idx.N_p + 3;
                DD = zeros(n, n);
                
                % Mesh stiffness:
                % Mesh component:
                sun_pla = stage_idx.sub_set('sun_planet');
                pla_rng = stage_idx.sub_set('planet_ring');

                d_sp = mesh_damp; %*sun_pla.k_mesh;
                d_rp = mesh_damp; %*pla_rng.k_mesh;
                
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
        
        function KK = stage_faulty_stiffness_matrix(stage_idx, val)

            if(strcmp(stage_idx.configuration, 'parallel'))
                n = 3;
                KK = zeros(n, n);
                
                r_p = (stage_idx.d(1)*1.0e-3)/2.0;
                r_w = (stage_idx.d(2)*1.0e-3)/2.0;
                
                k_pw = stage_idx.k_mesh;
                k_pw = k_pw*(1.0 - val);
                
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
                
                vec = ones(1, stage_idx.N_p);
                vec(2) = 1.0 - val;
                
                for jdx = 2:(n - 2)
                    
                    k_rp_idx = k_rp*(1.0 - vec(jdx - 1));
                    k_sp_idx = k_sp*(1.0 - vec(jdx - 1));
                    
                    KK(    1,   jdx) = r_c*r_p*(k_rp_idx - k_sp_idx);
                    KK(jdx  ,     1) = KK(1, jdx);
                    KK(jdx  ,   jdx) = (k_rp_idx + k_sp_idx)*r_p^2;
                    KK(n - 1,   jdx) = r_s*r_p*k_sp_idx;
                    KK(jdx  , n - 1) = KK(n - 1, jdx);
                end
            end
            range = (n - 1):n;
            
            KK(range, range) = KK(range, range) + ...
                stage_idx.output_shaft.stiffness_matrix('torsional');
            
        end
        
        function flag = test()
            g_set = NREL_5MW.gear_stage(1);
            MM       = Kahraman_94.stage_inertia_matrix(g_set);
            [Kb] = Kahraman_94.stage_stiffness_matrix(g_set);
            
            import matlab.mock.TestCase;
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
            dyn_for_mock.K_m = Kb*0.0;
            dyn_for_mock.M = MM;
            dyn_for_mock.A = [zeros(length(MM)),   eye(length(MM));
                                       -MM\Kb  , zeros(length(MM))]; 
            
            fn_test = dyn_for_mock.modal_analysis();
            fn_test = sort(fn_test);
            
            m_s = g_set.mass(1);    m_p = g_set.mass(2);    m_c = g_set.carrier.mass;
            n = g_set.N_p;
            k_1 = g_set.sub_set('planet_ring').k_mesh;
            k_2 = g_set.sub_set('sun_planet' ).k_mesh;
            
            % Fixed ring:
            Lambda_1 = m_p*m_c*m_s;
            Lambda_2 = -(n*k_2*m_p*m_c + (k_1 + k_2)*m_c*m_s + n*(k_1 + k_2)*m_p*m_s);
            Lambda_3 =   n*k_1*k_2*(n*m_p + m_c + 4.0*m_s);
            
            eig_val_1 = (-Lambda_2 - sqrt(Lambda_2^2 - 4.0*Lambda_1*Lambda_3))/(2.0*Lambda_1);
            eig_val_2 = (-Lambda_2 + sqrt(Lambda_2^2 - 4.0*Lambda_1*Lambda_3))/(2.0*Lambda_1);
            
            omega = [0.0;
                     repmat(sqrt((k_1 + k_2)/m_p), n - 1, 1);
                     sqrt(eig_val_1);
                     sqrt(eig_val_2)];
            fn_ana = omega./(2.0*pi);
            
        end
        
    end
end

