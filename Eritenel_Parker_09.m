classdef Eritenel_Parker_09 < Dynamic_Formulation
    %ERITENEL_PARKER_09 Calculates the inertia and stiffness matrices of a
    % multi-stage Drivetrain object according to:
    % [1] T. Eritenel and R. G. Parker, 'Modal properties of three-dimensional
    % helical planetary gears', J. Sound Vib., vol. 325, no. 1-2,
    % pp. 397-420, Aug. 2009. https://doi.org/10.1016/j.jsv.2009.03.002
    % [2] T. Eritenel, 'Three-Dimensional Nonlinear Dynamics and Vibration
    % Reduction of Gear Pairs and Planetary Gears', Ph.D Thesis,
    % Ohio State University, 2011.
    % http://rave.ohiolink.edu/etdc/view?acc_num=osu1298651902
    % [3] T. Eritenel and R. G. Parker, 'An investigation of tooth mesh
    % nonlinearity and partial contact loss in gear pairs using a
    % lumped-parameter model', Mech. Mach. Theory, vol. 56, pp. 28-51,
    % Oct. 2012. https://doi.org/10.1016/j.mechmachtheory.2012.05.002
    %
    % Rotation around z-axis.

    methods
        function obj = Eritenel_Parker_09(DT)
            obj@Dynamic_Formulation(DT);

            [obj.M, obj.G] = obj.inertia_matrix();
            obj.K = obj.stiffness_matrix();

        end
    end

    methods(Static)
        function [MM, GG] = stage_inertia_matrix(stage_idx)
            %
            
            G_sub = @(x)([0.0, -  x, 0.0;
                            x,  0.0, 0.0;
                          0.0,  0.0, 0.0]);
            
            if(strcmp(stage_idx.configuration, 'parallel'))
                % pinion, wheel
                % Inertia matrix:
                m_p = stage_idx.mass(1);
                m_w = stage_idx.mass(2);
                
                J_px = stage_idx.J_y(1);
                J_wx = stage_idx.J_y(2);
                
                J_pz = stage_idx.J_x(1);
                J_wz = stage_idx.J_x(2);
                
                diag_inertia = [J_px, J_px, J_pz, ones(1, 3)*m_p, ...
                                J_wx, J_wx, J_wz, ones(1, 3)*m_w];
                                
                % Gyroscopic matrix:
                r_p = (stage_idx.d_b(1)*1.0e-3)/2.0;
                r_w = (stage_idx.d_b(2)*1.0e-3)/2.0;

                G_p_rot = G_sub(r_p*(2.0*J_px - J_pz));
                G_p_tra = G_sub(2.0*m_p);
                G_w_rot = G_sub(r_w*(2.0*J_wx - J_wz));
                G_w_tra = G_sub(2.0*m_w);

                GG = blkdiag(G_p_rot, G_p_tra, ...
                             G_w_rot, G_w_tra);
                
            elseif(strcmp(stage_idx.configuration, 'planetary'))
                % sun, ring, carrier, N planets
                
                % Inertia matrix:
                m_s = stage_idx.mass(1);
                m_r = stage_idx.mass(3);
                m_c = stage_idx.carrier.mass;
                m_p = stage_idx.mass(2);
                
                J_sx = stage_idx.J_y(1);
                J_rx = stage_idx.J_y(3);
                J_cx = stage_idx.carrier.J_y;
                J_px = stage_idx.J_y(2);
                
                J_sz = stage_idx.J_x(1);
                J_rz = stage_idx.J_x(3);
                J_cz = stage_idx.carrier.J_x;
                J_pz = stage_idx.J_x(2);
                
                diag_inertia = [J_sx, J_sx, J_sz, ones(1, 3)*m_s, ...
                                J_rx, J_rx, J_rz, ones(1, 3)*m_r, ...
                                J_cx, J_cx, J_cz, ones(1, 3)*m_c, ...
                               [J_px, J_px, J_pz, ones(1, 3)*m_p]*repmat(eye(6), 1, stage_idx.N_p)];
                
                % Gyroscopic matrix:
                r_s = (stage_idx.d_b(1)*1.0e-3)/2.0;
                r_r = (stage_idx.d_b(3)*1.0e-3)/2.0;
                r_c =  stage_idx.a_w   *1.0e-3;
                r_p = (stage_idx.d_b(2)*1.0e-3)/2.0;
                
                G_s_rot = G_sub(r_s*(2.0*J_sx - J_sz));
                G_s_tra = G_sub(2.0*m_s);
                G_r_rot = G_sub(r_r*(2.0*J_rx - J_rz));
                G_r_tra = G_sub(2.0*m_r);
                G_c_rot = G_sub(r_c*(2.0*J_cx - J_cz));
                G_c_tra = G_sub(2.0*m_c);
                G_p_rot = G_sub(r_p*(2.0*J_px - J_pz));
                G_p_tra = G_sub(2.0*m_p);
                
                G_p = blkdiag(G_p_rot, G_p_tra);
                G_p = repmat({G_p}, 1, stage_idx.N_p);
                G_p = blkdiag(G_p{:});
                
                GG = blkdiag(G_s_rot, G_s_tra, ...
                             G_r_rot, G_r_tra, ...
                             G_c_rot, G_c_tra, ...
                             G_p);
            end
            
            MM = diag(diag_inertia);
        end

        function KK = stage_stiffness_matrix(stage_idx)

            % Eq. (29):
            R = [0, 1, 0, 0, 0, 0; ...
                -1, 0, 0, 0, 0, 0; ...
                0, 0, 1, 0, 0, 0; ...
                0, 0, 0, 0, 1, 0; ...
                0, 0, 0, -1, 0, 0; ...
                0, 0, 0, 0, 0, 1];

            if (strcmp(stage_idx.configuration, 'parallel'))
            elseif (strcmp(stage_idx.configuration, 'planetary'))
                N1 = stage_idx.N_p - 1;

                sun_pla = stage_idx.sub_set('sun_planet');
                pla_rng = stage_idx.sub_set('planet_ring');

                % OBS.: got a negative k_mesh for DTU_10MW planet-ring
                k_sp = sun_pla.k_mesh;
                k_rp = pla_rng.k_mesh;
                
                % Translational gear mesh stiffness:
                % - odd index: sun-planet;
                % - even index: ring-planet;
                k = repmat([k_sp, k_rp], 1, stage_idx.N_p);

                % Angular position of planets:
                alpha = (0:N1)*(360.0/stage_idx.N_p);
                % Base helix angle:
                psi = stage_idx.beta_b;
                % Transverse operating pressure angle:
                phi_sp = sun_pla.alpha_wt;
                phi_rp = pla_rng.alpha_wt;

                sum_sin = sum(sind(alpha));
                sum_cos = sum(cosd(alpha));
                sum_sin2 = sum(sind(alpha).^2);
                sum_cos2 = sum(cosd(alpha).^2);
                sum_sc = sum(sind(alpha).*cosd(alpha));

                r_s = (stage_idx.d_b(1)*1.0e-3)/2.0;
                r_r = (stage_idx.d_b(3)*1.0e-3)/2.0;
%                 r_c =  stage_idx.a_w   *1.0e-3;
                r_p = (stage_idx.d_b(2)*1.0e-3)/2.0;
                
                % Center of stiffness (axial position):
                c = zeros(size(k));
                
                % Center of mass (axial position):
                e_s = 0.0;
                e_r = e_s;
                e_c = e_s;
                e_p = e_s;
                
                % tilting stiffness of bearing and mesh:
                kappa = k*0;
                kappa_As  = stage_idx.bearing(5).K_beta;
                kappa_Azs = stage_idx.bearing(5).K_alpha;
                kappa_Bs  = kappa_As;
                kappa_Bzs = kappa_Azs;
                
                kappa_Ar = 0.0;
                
                % Bearing axial position:
                L_As = 0.0;                L_Bs = 0.0;
                L_Ac = 0.0;                L_Bc = 0.0;
                L_Ap = 0.0;                L_Bp = 0.0;
                L_Ar = 0.0;                L_Br = 0.0;

                % Bearing translational stiffness:
                k_As = stage_idx.bearing(5).K_y;        k_Bs = k_As;
                k_Ac = stage_idx.bearing(3).K_y;        k_Bc = stage_idx.bearing(4).K_y;
                k_Ap = stage_idx.bearing(1).K_y;        k_Bp = stage_idx.bearing(2).K_y;
                k_Ar = stage_idx.bearing(6).K_y;        k_Br = k_Ar;
                
                % Bearing axial stiffness:
                k_Azs = stage_idx.bearing(5).K_x;
                k_Bzs = k_Azs;

                %% Sun gear matrices:
                % Eq. (87):
                Upsi_s = @(jdx)([0 0 k(jdx) * r_s * ((e_s - c(jdx)) * cosd(psi) - r_s * tand(phi_sp) * sind(psi)) * cosd(psi) 0 0 -k(jdx) * ((e_s - c(jdx)) * cosd(psi) - r_s * tand(phi_sp) * sind(psi)) * sind(psi); ...
                    0 0 k(jdx) * r_s ^ 2 * sind(psi) * cosd(psi) 0 0 -k(jdx) * r_s * sind(psi) ^ 2; ...
                    0 0 0 0 k(jdx) * r_s * cosd(psi) ^ 2 0; ...
                    0 0 0 0 0 0; ...
                    0 0 0 0 0 0; ...
                    0 0 0 0 0 0]);

                % Eq. (88):
                Theta_s = @(jdx)([kappa(jdx) + k(jdx) * ((e_s - c(jdx)) * cosd(psi) - r_s * tand(phi_sp) * sind(psi)) ^ 2 k(jdx) * ((e_s - c(jdx)) * cosd(psi) - r_s * tand(phi_sp) * sind(psi)) * r_s * sind(psi) 0 0 k(jdx) * ((e_s - c(jdx)) * cosd(psi) - r_s * tand(phi_sp) * sind(psi)) * cosd(psi) 0; ...
                    0 k(jdx) * r_s * sind(psi) ^ 2 0 0 k(jdx) * r_s * sind(psi) * cosd(psi) 0; ...
                    0 0 0 0 0 0; ...
                    0 0 0 0 0 0; ...
                    0 0 0 0 k(jdx) * cosd(psi) ^ 2 0; ...
                    0 0 0 0 0 0]);

                % Eq. (89):
                Xi_s = @(jdx)([kappa(jdx) + k(jdx) * ((e_s - c(jdx)) * cosd(psi) - r_s * tand(phi_sp) * sind(psi)) ^ 2 k(jdx) * ((e_s - c(jdx)) * cosd(psi) - r_s * tand(phi_sp) * sind(psi)) * r_s * sind(psi) 0 0 k(jdx) * ((e_s - c(jdx)) * cosd(psi) - r_s * tand(phi_sp) * sind(psi)) * cosd(psi) 0; ...
                    0 k(jdx) * r_s * sind(psi) ^ 2 0 0 k(jdx) * r_s * sind(psi) * cosd(psi) 0; ...
                    0 0 0 0 0 0; ...
                    0 0 0 0 0 0; ...
                    0 0 0 0 k(jdx) * cosd(psi) ^ 2 0; ...
                    0 0 0 0 0 0]);

                % Eq. (90):
                Psi_s = @(jdx)([(-L_As - e_s) ^ 2 * k_As + (L_Bs - e_s) ^ 2 * k_Bs + kappa_As + kappa_Bs 0 0 0 -(-L_As - e_s) * k_As - (L_Bs - e_s) * k_Bs 0; ...
                    0 (-L_As - e_s) ^ 2 * k_As + (L_Bs - e_s) ^ 2 * k_Bs + kappa_As + kappa_Bs 0 (-L_As - e_s) * k_As + (L_Bs - e_s) * k_Bs 0 0; ...
                    0 0 kappa_Azs + kappa_Bzs + k(jdx) * r_s ^ 2 * cosd(psi) ^ 2 0 0 -k(jdx) * r_s * sind(psi) * cosd(psi); ...
                    0 0 0 k_As + k_Bs 0 0; ...
                    0 0 0 0 k_As + k_Bs 0; ...
                    0 0 0 0 0 k_Azs + k_Bzs + k(jdx) * sind(psi) ^ 2]);

                % Eq. (91):
                Lambda_s = @(jdx)([k(jdx) * ((e_s - c(jdx)) * cosd(psi) - r_s * tand(phi_sp)* sind(psi)) * r_p * sind(psi) k(jdx) * ((e_s - c(jdx)) * cosd(psi) - r_s * tand(phi_sp)* sind(psi)) * ((e_p - c(jdx)) * cosd(psi) + r_p * tand(phi_sp)* sind(psi)) + kappa(jdx) k(jdx) * ((e_s - c(jdx)) * cosd(psi) - r_s * tand(phi_sp)* sind(psi)) * r_p * cosd(psi) -k(jdx) * ((e_s - c(jdx)) * cosd(psi) - r_s * tand(phi_sp)* sind(psi)) * cosd(psi) 0 k(jdx) * ((e_s - c(jdx)) * cosd(psi) - r_s * tand(phi_sp)* sind(psi)) * sind(psi); ...
                    k(jdx) * r_s * r_p * sind(psi) ^ 2 k(jdx) * r_s * sind(psi) * ((e_p - c(jdx)) * cosd(psi) + r_p * tand(phi_sp)* sind(psi)) k(jdx) * r_s * r_p * sind(psi) * cosd(psi) -k(jdx) * r_s * sind(psi) * cosd(psi) 0 k(jdx) * r_s * sind(psi) ^ 2; ...
                    0 0 0 0 0 0; ...
                    0 0 0 0 0 0; ...
                    k(jdx) * r_p * sind(psi) * cosd(psi) k(jdx) * ((e_p - c(jdx)) * cosd(psi) + r_p * tand(phi_sp)* sind(psi)) * cosd(psi) k(jdx) * r_p * cosd(psi) ^ 2 -k(jdx) * cosd(psi) ^ 2 0 k(jdx) * sind(psi) * cosd(psi); ...
                    0 0 0 0 0 0]);

                % Eq. (92):
                Gamma_s = @(jdx)([0 0 0 0 0 0; ...
                    0 0 0 0 0 0; ...
                    k(jdx) * r_s * r_p * sind(psi) * cosd(psi) k(jdx) * r_s * cosd(psi) * ((e_p - c(jdx)) * cosd(psi) + r_p * tand(phi_sp)* sind(psi)) k(jdx) * r_s * r_p * cosd(psi) ^ 2 -k(jdx) * r_s * cosd(psi) ^ 2 0 k(jdx) * r_s * sind(psi) * cosd(psi); ...
                    0 0 0 0 0 0; ...
                    0 0 0 0 0 0; ...
                    -k(jdx) * r_p * sind(psi) ^ 2 -k(jdx) * sind(psi) * ((e_p - c(jdx)) * cosd(psi) + r_p * tand(phi_sp)* sind(psi)) -k(jdx) * r_p * sind(psi) * cosd(psi) k(jdx) * sind(psi) * cosd(psi) 0 -k(jdx) * sind(psi) ^ 2]);

                % Eq. (28):
                K_s = Upsi_s(1) * sum_sin + R * Upsi_s(1) * R' * sum_cos + Theta_s(1) * sum_sin2 + ...
                    R * Theta_s(1) * R' * sum_cos2 + Xi_s(1) * sum_sc + Psi_s(1);

                % Eq. (31):
                K_si = @(i)(Lambda_s(i) * sind(alpha(i)) + R * Lambda_s(i) * cosd(alpha(i)) + Gamma_s(i));

                %% Ring gear matrices:
                % Eq. (93):
                Upsi_r = @(jdx)([0 0 k(jdx) * r_r * (cosd(phi_sp + phi_rp) * ((e_r - c(jdx)) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * sind(phi_sp + phi_rp)) * cosd(psi) 0 0 k(jdx) * (cosd(phi_sp + phi_rp) * ((e_r - c(jdx)) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * sind(phi_sp + phi_rp)) * sind(psi);
                    0 0 k(jdx) * r_r * (sind(phi_sp + phi_rp) * ((e_r - c(jdx)) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * cosd(phi_sp + phi_rp)) * cosd(psi) 0 0 k(jdx) * (sind(phi_sp + phi_rp) * ((e_r - c(jdx)) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * cosd(phi_sp + phi_rp)) * sind(psi);
                    0 0 0 -k(jdx) * r_r * cosd(psi) ^ 2 * sind(phi_sp + phi_rp) k(jdx) * r_r * cosd(psi) ^ 2 * cosd(phi_sp + phi_rp) 0;
                    0 0 0 0 0 -k(jdx) * cosd(psi) * sind(phi_sp + phi_rp) * sind(psi);
                    0 0 0 0 0 k(jdx) * cosd(psi) * cosd(phi_sp + phi_rp) * sind(psi);
                    0 0 0 0 0 0]);
                
                % Eq. (94):
                Theta_r = @(jdx)([kappa(jdx) * cosd(phi_sp + phi_rp) ^ 2 + k(jdx) * (cosd(phi_sp + phi_rp) * ((e_r - c(jdx)) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * sind(phi_sp + phi_rp)) ^ 2 kappa(jdx) * cosd(phi_sp + phi_rp) * sind(phi_sp + phi_rp) + k(jdx) * (cosd(phi_sp + phi_rp) * ((e_r - c(jdx)) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * sind(phi_sp + phi_rp)) * (sind(phi_sp + phi_rp) * ((e_r - c(jdx)) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * cosd(phi_sp + phi_rp)) 0 -k(jdx) * (cosd(phi_sp + phi_rp) * ((e_r - c(jdx)) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * sind(phi_sp + phi_rp)) * cosd(psi) * sind(phi_sp + phi_rp) k(jdx) * (cosd(phi_sp + phi_rp) * ((e_r - c(jdx)) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * sind(phi_sp + phi_rp)) * cosd(psi) * cosd(phi_sp + phi_rp) 0;
                    0 kappa(jdx) * sind(phi_sp + phi_rp) ^ 2 + k(jdx) * (sind(phi_sp + phi_rp) * ((e_r - c(jdx)) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * cosd(phi_sp + phi_rp)) ^ 2 0 -k(jdx) * (sind(phi_sp + phi_rp) * ((e_r - c(jdx)) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * cosd(phi_sp + phi_rp)) * cosd(psi) * sind(phi_sp + phi_rp) k(jdx) * (sind(phi_sp + phi_rp) * ((e_r - c(jdx)) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * cosd(phi_sp + phi_rp)) * cosd(psi) * cosd(phi_sp + phi_rp) 0;
                    0 0 0 0 0 0;
                    0 0 0 k(jdx) * cosd(psi) ^ 2 * sind(phi_sp + phi_rp) ^ 2 -k(jdx) * cosd(psi) ^ 2 * sind(phi_sp + phi_rp) * cosd(phi_sp + phi_rp) 0;
                    0 0 0 0 k(jdx) * cosd(psi) ^ 2 * cosd(phi_sp + phi_rp) ^ 2 0;
                    0 0 0 0 0 0]);
                
                % Eq. (95):
                Xi_r = @(jdx)([0.2e1 * kappa(jdx) * cosd(phi_sp + phi_rp) * sind(phi_sp + phi_rp) + 0.2e1 * k(jdx) * (cosd(phi_sp + phi_rp) * ((e_r - c(jdx)) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * sind(phi_sp + phi_rp)) * (sind(phi_sp + phi_rp) * ((e_r - c(jdx)) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * cosd(phi_sp + phi_rp)) k(jdx) * ((sind(phi_sp + phi_rp) * ((e_r - c(jdx)) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * cosd(phi_sp + phi_rp)) ^ 2 - (cosd(phi_sp + phi_rp) * ((e_r - c(jdx)) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * sind(phi_sp + phi_rp)) ^ 2) + kappa(jdx) * (sind(phi_sp + phi_rp) ^ 2 - cosd(phi_sp + phi_rp) ^ 2) 0 k(jdx) * ((cosd(phi_sp + phi_rp) * ((e_r - c(jdx)) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * sind(phi_sp + phi_rp)) * cosd(psi) * cosd(phi_sp + phi_rp) - (sind(phi_sp + phi_rp) * ((e_r - c(jdx)) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * cosd(phi_sp + phi_rp)) * cosd(psi) * sind(phi_sp + phi_rp)) k(jdx) * ((sind(phi_sp + phi_rp) * ((e_r - c(jdx)) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * cosd(phi_sp + phi_rp)) * cosd(psi) * cosd(phi_sp + phi_rp) - (cosd(phi_sp + phi_rp) * ((e_r - c(jdx)) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * sind(phi_sp + phi_rp)) * cosd(psi) * sind(phi_sp + phi_rp)) 0;
                    0 -0.2e1 * kappa(jdx) * cosd(phi_sp + phi_rp) * sind(phi_sp + phi_rp) - 0.2e1 * k(jdx) * (cosd(phi_sp + phi_rp) * ((e_r - c(jdx)) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * sind(phi_sp + phi_rp)) * (sind(phi_sp + phi_rp) * ((e_r - c(jdx)) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * cosd(phi_sp + phi_rp)) 0 k(jdx) * ((sind(phi_sp + phi_rp) * ((e_r - c(jdx)) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * cosd(phi_sp + phi_rp)) * cosd(psi) * cosd(phi_sp + phi_rp) - (cosd(phi_sp + phi_rp) * ((e_r - c(jdx)) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * sind(phi_sp + phi_rp)) * cosd(psi) * sind(phi_sp + phi_rp)) -k(jdx) * ((cosd(phi_sp + phi_rp) * ((e_r - c(jdx)) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * sind(phi_sp + phi_rp)) * cosd(psi) * cosd(phi_sp + phi_rp) - (sind(phi_sp + phi_rp) * ((e_r - c(jdx)) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * cosd(phi_sp + phi_rp)) * cosd(psi) * sind(phi_sp + phi_rp)) 0;
                    0 0 0 0 0 0;
                    0 0 0 -0.2e1 * k(jdx) * cosd(psi) ^ 2 * sind(phi_sp + phi_rp) * cosd(phi_sp + phi_rp) k(jdx) * (-cosd(psi) ^ 2 * sind(phi_sp + phi_rp) ^ 2 + cosd(psi) ^ 2 * cosd(phi_sp + phi_rp) ^ 2) 0;
                    0 0 0 0 0.2e1 * k(jdx) * cosd(psi) ^ 2 * sind(phi_sp + phi_rp) * cosd(phi_sp + phi_rp) 0;
                    0 0 0 0 0 0]);
                
                % Eq. (96):
                Psi_r = @(jdx)([k_Ar * (-L_Ar - e_r) ^ 2 + k_Br * (L_Br - e_r) ^ 2 + kappa_Ar + kappa_Br  0   0   0   -(-L_Ar - e_r) * k_Ar - (L_Br - e_r) * k_Br   0;
                    0   k_Ar * (-L_Ar - e_r) ^ 2 + k_Br * (L_Br - e_r) ^ 2 + kappa_Ar + kappa_Br  0   (-L_Ar - e_r) * k_Ar + (L_Br - e_r) * k_Br    0   0;
                    0   0   kappa_Azr + kappa_Bzr + k(jdx) * RR ^ 2 * cosd(psi) ^ 2  0   0   k(jdx) * r_r * cosd(psi) * sind(psi);
                    0   (-L_Ar - e_r) * k_Ar + (L_Br - e_r) * k_Br    0   k_Ar + k_Br 0   0;
                    -(-L_Ar - e_r) * k_Ar - (L_Br - e_r) * k_Br   0   0   0   k_Ar + k_Br 0;
                    0   0   k(jdx) * r_r * cosd(psi) * sind(psi) 0   0   k_Azr + k_Bzr + k(jdx) * sind(psi) ^ 2]);
                
                % Eq. (97):
                Lambda_r = @(jdx)([k(jdx) * (cosd(phi_sp + phi_rp) * ((e_r - c(jdx)) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * sind(phi_sp + phi_rp)) * (sind(phi_sp + phi_rp) * ((e_p - c(jdx)) * cosd(psi) - r_p * tand(phi_rp) * sind(psi)) + r_p * sind(psi) * cosd(phi_sp + phi_rp)) - kappa(jdx) * cosd(phi_sp + phi_rp) * sind(phi_sp + phi_rp) k(jdx) * (cosd(phi_sp + phi_rp) * ((e_r - c(jdx)) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * sind(phi_sp + phi_rp)) * (cosd(phi_sp + phi_rp) * ((e_p - c(jdx)) * cosd(psi) - r_p * tand(phi_rp) * sind(psi)) + r_p * sind(psi) * sind(phi_sp + phi_rp)) - kappa(jdx) * cosd(phi_sp + phi_rp) ^ 2 -k(jdx) * r_p * (cosd(phi_sp + phi_rp) * ((e_r - c(jdx)) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * sind(phi_sp + phi_rp)) * cosd(psi) -k(jdx) * (cosd(phi_sp + phi_rp) * ((e_r - c(jdx)) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * sind(phi_sp + phi_rp)) * cosd(psi) * cosd(phi_sp + phi_rp) -k(jdx) * (cosd(phi_sp + phi_rp) * ((e_r - c(jdx)) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * sind(phi_sp + phi_rp)) * cosd(psi) * sind(phi_sp + phi_rp) -k(jdx) * (cosd(phi_sp + phi_rp) * ((e_r - c(jdx)) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * sind(phi_sp + phi_rp)) * sind(psi);
                    k(jdx) * (sind(phi_sp + phi_rp) * ((e_r - c(jdx)) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * cosd(phi_sp + phi_rp)) * (sind(phi_sp + phi_rp) * ((e_p - c(jdx)) * cosd(psi) - r_p * tand(phi_rp) * sind(psi)) + r_p * sind(psi) * cosd(phi_sp + phi_rp)) - kappa(jdx) * sind(phi_sp + phi_rp) ^ 2 k(jdx) * (sind(phi_sp + phi_rp) * ((e_r - c(jdx)) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * cosd(phi_sp + phi_rp)) * (cosd(phi_sp + phi_rp) * ((e_p - c(jdx)) * cosd(psi) - r_p * tand(phi_rp) * sind(psi)) + r_p * sind(psi) * sind(phi_sp + phi_rp)) + kappa(jdx) * cosd(phi_sp + phi_rp) * sind(phi_sp + phi_rp) -k(jdx) * (sind(phi_sp + phi_rp) * ((e_r - c(jdx)) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * cosd(phi_sp + phi_rp)) * r_p * cosd(psi) -k(jdx) * (sind(phi_sp + phi_rp) * ((e_r - c(jdx)) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * cosd(phi_sp + phi_rp)) * cosd(psi) * cosd(phi_sp + phi_rp) -k(jdx) * (sind(phi_sp + phi_rp) * ((e_r - c(jdx)) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * cosd(phi_sp + phi_rp)) * cosd(psi) * sind(phi_sp + phi_rp) -k(jdx) * (sind(phi_sp + phi_rp) * ((e_r - c(jdx)) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * cosd(phi_sp + phi_rp)) * sind(psi);
                    0 0 0 0 0 0;
                    -k(jdx) * (sind(phi_sp + phi_rp) * ((e_p - c(jdx)) * cosd(psi) - r_p * tand(phi_rp) * sind(psi)) + r_p * sind(psi) * cosd(phi_sp + phi_rp)) * cosd(psi) * sind(phi_sp + phi_rp) -k(jdx) * (cosd(phi_sp + phi_rp) * ((e_p - c(jdx)) * cosd(psi) - r_p * tand(phi_rp) * sind(psi)) + r_p * sind(psi) * sind(phi_sp + phi_rp)) * cosd(psi) * sind(phi_sp + phi_rp) -k(jdx) * cosd(psi) ^ 2 * sind(phi_sp + phi_rp) * r_p k(jdx) * cosd(psi) ^ 2 * sind(phi_sp + phi_rp) * cosd(phi_sp + phi_rp) k(jdx) * cosd(psi) ^ 2 * sind(phi_sp + phi_rp) ^ 2 k(jdx) * cosd(psi) * sind(phi_sp + phi_rp) * sind(psi);
                    k(jdx) * (sind(phi_sp + phi_rp) * ((e_p - c(jdx)) * cosd(psi) - r_p * tand(phi_rp) * sind(psi)) + r_p * sind(psi) * cosd(phi_sp + phi_rp)) * cosd(psi) * cosd(phi_sp + phi_rp) k(jdx) * (cosd(phi_sp + phi_rp) * ((e_p - c(jdx)) * cosd(psi) - r_p * tand(phi_rp) * sind(psi)) + r_p * sind(psi) * sind(phi_sp + phi_rp)) * cosd(psi) * cosd(phi_sp + phi_rp) k(jdx) * cosd(psi) ^ 2 * cosd(phi_sp + phi_rp) * r_p -k(jdx) * cosd(psi) ^ 2 * cosd(phi_sp + phi_rp) ^ 2 -k(jdx) * cosd(psi) ^ 2 * sind(phi_sp + phi_rp) * cosd(phi_sp + phi_rp) -k(jdx) * cosd(psi) * cosd(phi_sp + phi_rp) * sind(psi);
                    0 0 0 0 0 0]);
                
                % Eq. (98):
                Gamma_r = @(jdx)([k(jdx) * (cosd(phi_sp + phi_rp) * ((e_r - c(jdx)) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * sind(phi_sp + phi_rp)) * (sind(phi_sp + phi_rp) * ((e_p - c(jdx)) * cosd(psi) - r_p * tand(phi_rp) * sind(psi)) + r_p * sind(psi) * cosd(phi_sp + phi_rp)) - kappa(jdx) * cosd(phi_sp + phi_rp) * sind(phi_sp + phi_rp) k(jdx) * (cosd(phi_sp + phi_rp) * ((e_r - c(jdx)) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * sind(phi_sp + phi_rp)) * (cosd(phi_sp + phi_rp) * ((e_p - c(jdx)) * cosd(psi) - r_p * tand(phi_rp) * sind(psi)) + r_p * sind(psi) * sind(phi_sp + phi_rp)) - kappa(jdx) * cosd(phi_sp + phi_rp) ^ 2 k(jdx) * r_p * (cosd(phi_sp + phi_rp) * ((e_r - c(jdx)) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * sind(phi_sp + phi_rp)) * cosd(psi) -k(jdx) * (cosd(phi_sp + phi_rp) * ((e_r - c(jdx)) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * sind(phi_sp + phi_rp)) * cosd(psi) * cosd(phi_sp + phi_rp) -k(jdx) * (cosd(phi_sp + phi_rp) * ((e_r - c(jdx)) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * sind(phi_sp + phi_rp)) * cosd(psi) * sind(phi_sp + phi_rp) -k(jdx) * (cosd(phi_sp + phi_rp) * ((e_r - c(jdx)) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * sind(phi_sp + phi_rp)) * sind(psi);
                    k(jdx) * (sind(phi_sp + phi_rp) * ((e_r - c(jdx)) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * cosd(phi_sp + phi_rp)) * (sind(phi_sp + phi_rp) * ((e_p - c(jdx)) * cosd(psi) - r_p * tand(phi_rp) * sind(psi)) + r_p * sind(psi) * cosd(phi_sp + phi_rp)) - kappa(jdx) * sind(phi_sp + phi_rp) ^ 2 k(jdx) * (sind(phi_sp + phi_rp) * ((e_r - c(jdx)) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * cosd(phi_sp + phi_rp)) * (cosd(phi_sp + phi_rp) * ((e_p - c(jdx)) * cosd(psi) - r_p * tand(phi_rp) * sind(psi)) + r_p * sind(psi) * sind(phi_sp + phi_rp)) + kappa(jdx) * cosd(phi_sp + phi_rp) * sind(phi_sp + phi_rp) 0 -k(jdx) * (sind(phi_sp + phi_rp) * ((e_r - c(jdx)) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * cosd(phi_sp + phi_rp)) * cosd(psi) * cosd(phi_sp + phi_rp) -k(jdx) * (sind(phi_sp + phi_rp) * ((e_r - c(jdx)) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * cosd(phi_sp + phi_rp)) * cosd(psi) * sind(phi_sp + phi_rp) -k(jdx) * (sind(phi_sp + phi_rp) * ((e_r - c(jdx)) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * cosd(phi_sp + phi_rp)) * sind(psi);
                    k(jdx) * r_r * cosd(psi) * (sind(phi_sp + phi_rp) * ((e_p - c(jdx)) * cosd(psi) - r_p * tand(phi_rp) * sind(psi)) + r_p * sind(psi) * cosd(phi_sp + phi_rp)) k(jdx) * r_r * cosd(psi) * (cosd(phi_sp + phi_rp) * ((e_p - c(jdx)) * cosd(psi) - r_p * tand(phi_rp) * sind(psi)) + r_p * sind(psi) * sind(phi_sp + phi_rp)) -k(jdx) * r_r * r_p * cosd(psi) ^ 2 -k(jdx) * r_r * cosd(psi) ^ 2 * cosd(phi_sp + phi_rp) -k(jdx) * r_r * cosd(psi) ^ 2 * sind(phi_sp + phi_rp) -k(jdx) * r_r * cosd(psi) * sind(psi);
                    -k(jdx) * (sind(phi_sp + phi_rp) * ((e_p - c(jdx)) * cosd(psi) - r_p * tand(phi_rp) * sind(psi)) + r_p * sind(psi) * cosd(phi_sp + phi_rp)) * cosd(psi) * sind(phi_sp + phi_rp) -k(jdx) * (cosd(phi_sp + phi_rp) * ((e_p - c(jdx)) * cosd(psi) - r_p * tand(phi_rp) * sind(psi)) + r_p * sind(psi) * sind(phi_sp + phi_rp)) * cosd(psi) * sind(phi_sp + phi_rp) k(jdx) * cosd(psi) ^ 2 * sind(phi_sp + phi_rp) * r_p k(jdx) * cosd(psi) ^ 2 * sind(phi_sp + phi_rp) * cosd(phi_sp + phi_rp) k(jdx) * cosd(psi) ^ 2 * sind(phi_sp + phi_rp) ^ 2 k(jdx) * cosd(psi) * sind(phi_sp + phi_rp) * sind(psi);
                    k(jdx) * (sind(phi_sp + phi_rp) * ((e_p - c(jdx)) * cosd(psi) - r_p * tand(phi_rp) * sind(psi)) + r_p * sind(psi) * cosd(phi_sp + phi_rp)) * cosd(psi) * cosd(phi_sp + phi_rp) k(jdx) * (cosd(phi_sp + phi_rp) * ((e_p - c(jdx)) * cosd(psi) - r_p * tand(phi_rp) * sind(psi)) + r_p * sind(psi) * sind(phi_sp + phi_rp)) * cosd(psi) * cosd(phi_sp + phi_rp) k(jdx) * cosd(psi) ^ 2 * cosd(phi_sp + phi_rp) * r_p -k(jdx) * cosd(psi) ^ 2 * cosd(phi_sp + phi_rp) ^ 2 -k(jdx) * cosd(psi) ^ 2 * sind(phi_sp + phi_rp) * cosd(phi_sp + phi_rp) -k(jdx) * cosd(psi) * cosd(phi_sp + phi_rp) * sind(psi);
                    k(jdx) * sind(psi) * (sind(phi_sp + phi_rp) * ((e_p - c(jdx)) * cosd(psi) - r_p * tand(phi_rp) * sind(psi)) + r_p * sind(psi) * cosd(phi_sp + phi_rp)) k(jdx) * sind(psi) * (cosd(phi_sp + phi_rp) * ((e_p - c(jdx)) * cosd(psi) - r_p * tand(phi_rp) * sind(psi)) + r_p * sind(psi) * sind(phi_sp + phi_rp)) -k(jdx) * r_p * sind(psi) * cosd(psi) -k(jdx) * cosd(psi) * cosd(phi_sp + phi_rp) * sind(psi) -k(jdx) * cosd(psi) * sind(phi_sp + phi_rp) * sind(psi) -k(jdx) * r_r * cosd(psi) * sind(psi)]);

                % Eq. (28):
                K_r = Upsi_r(2) * sum_sin + R * Upsi_r(2) * R' * sum_cos + Theta_r(2) * sum_sin2 + ...
                    R * Theta_r(2) * R' * sum_cos2 + Xi_r(2) * sum_sc + Psi_r(2);

                % Eq. (31):
                K_ri = @(i)(Lambda_r(i) * sind(alpha(i)) + R * Lambda_r(i) * cosd(alpha(i)) + Gamma_r(i));

                %% Carrier matrices
                % Eq. (99):
                Upsi_c = [0, 0, -(r_s + r_p) * ((-L_Ap - e_s) * k_Ap + (L_Bp - e_s) * k_Bp), 0, 0, tand(phi_sp)* (r_s + r_p) * (k_Azp + k_Bzp); ...
                    0, 0, -tand(phi_sp)* (r_s + r_p) * ((-L_Ap - e_s) * k_Ap + (L_Bp - e_s) * k_Bp), 0, 0, -(r_s + r_p) * (k_Azp + k_Bzp); ...
                    0, 0, 0, -tand(phi_sp)* (r_s + r_p) * (k_Ap + k_Bp), (r_s + r_p) * (k_Ap + k_Bp), 0; ...
                    0, 0, 0, 0, 0, 0; ...
                    0, 0, 0, 0, 0, 0; ...
                    0, 0, 0, 0, 0, 0];

                % Eq. (100):
                Theta_c = [k_Ap * (-L_Ap - e_s) ^ 2 + k_Bp * (L_Bp - e_s) ^ 2 + (k_Azp + k_Bzp) * tand(phi_sp)^ 2 * (r_s + r_p) ^ 2 + kappa_Ap + kappa_Bp, -tand(phi_sp)* (r_s + r_p) ^ 2 * (k_Azp + k_Bzp), 0, 0, -(-L_Ap - e_s) * k_Ap - (L_Bp - e_s) * k_Bp, 0; ...
                    0, k_Ap * (-L_Ap - e_s) ^ 2 + k_Bp * (L_Bp - e_s) ^ 2 + (k_Azp + k_Bzp) * (r_s + r_p) ^ 2 + kappa_Ap + kappa_Bp, 0, (-L_Ap - e_s) * k_Ap + (L_Bp - e_s) * k_Bp, 0, 0; ...
                    0, 0, 0, 0, 0, 0; ...
                    0, 0, 0, k_Ap + k_Bp, 0, 0; ...
                    0, 0, 0, 0, k_Ap + k_Bp, 0; ...
                    0, 0, 0, 0, 0, 0];

                % Eq. (101):
                Xi_c = [-0.2e1 * tand(phi_sp)* (r_s + r_p) ^ 2 * (k_Azp + k_Bzp), (-tand(phi_sp)^2 * (r_s + r_p)^2 + (r_s + r_p)^2) * (k_Azp + k_Bzp), 0, 0, 0, 0; ...
                    0, 0.2e1 * tand(phi_sp)* (r_s + r_p) ^ 2 * (k_Azp + k_Bzp), 0, 0, 0, 0; ...
                    0, 0, 0, 0, 0, 0; ...
                    0, 0, 0, 0, 0, 0; ...
                    0, 0, 0, 0, 0, 0; ...
                    0, 0, 0, 0, 0, 0];

                % Eq. (102):
                Psi_c = [(-L_Ac - e_c) ^ 2 * k_Ac + (L_Bc - e_c) ^ 2 * k_Bc + kappa_Ac + kappa_Bc, 0, 0, 0, -(-L_Ac - e_c) * k_Ac - (L_Bc - e_c) * k_Bc, 0; ...
                    0, (-L_Ac - e_c) ^ 2 * k_Ac + (L_Bc - e_c) ^ 2 * k_Bc + kappa_Ac + kappa_Bc, 0, (-L_Ac - e_c) * k_Ac + (L_Bc - e_c) * k_Bc, 0, 0; ...
                    0, 0, kappa_Azc + kappa_Bzc + kappa_Azp + kappa_Bzp + (tand(phi_sp)^2 * (r_s + r_p)^2 + (r_s + r_p)^2) * (k_Ap + k_Bp), 0, 0, 0; ...
                    0, 0, 0, k_Ac + k_Bc, 0, 0; ...
                    0, 0, 0, 0, k_Ac + k_Bc, 0; ...
                    0, 0, 0, 0, 0, k_Azc + k_Bzc + k_Azp + k_Bzp];

                % Eq. (103):
                Lambda_c = [0, (-L_Ap - e_s) * (-L_Ap - e_p) * k_Ap + (L_Bp - e_s) * (L_Bp - e_p) * k_Bp + kappa_Ap + kappa_Bp, 0, (-L_Ap - e_s) * k_Ap + (L_Bp - e_s) * k_Bp, 0, -tand(phi_sp)* (r_s + r_p) * (k_Azp + k_Bzp); ...
                    -(-L_Ap - e_s) * (-L_Ap - e_p) * k_Ap - (L_Bp - e_s) * (L_Bp - e_p) * k_Bp - kappa_Ap - kappa_Bp, 0, 0, 0, (-L_Ap - e_s) * k_Ap + (L_Bp - e_s) * k_Bp, (r_s + r_p) * (k_Azp + k_Bzp); ...
                    tand(phi_sp)* (r_s + r_p) * ((-L_Ap - e_p) * k_Ap + (L_Bp - e_p) * k_Bp), -(r_s + r_p) * ((-L_Ap - e_p) * k_Ap + (L_Bp - e_p) * k_Bp), -kappa_Azp - kappa_Bzp, -(r_s + r_p) * (k_Ap + k_Bp), -tand(phi_sp)* (r_s + r_p) * (k_Ap + k_Bp), 0; ...
                    -(-L_Ap - e_p) * k_Ap - (L_Bp - e_p) * k_Bp, 0, 0, 0, k_Ap + k_Bp, 0; ...
                    0, -(-L_Ap - e_p) * k_Ap - (L_Bp - e_p) * k_Bp, 0, -k_Ap - k_Bp, 0, 0; ...
                    0, 0, 0, 0, 0, -k_Azp - k_Bzp];

                % Eq. (104):
                Gamma_c = [0 0 0 0 0 0; ...
                    0 0 0 0 0 0; ...
                    tand(phi_sp)* (r_s + r_p) * ((-L_Ap - e_p) * k_Ap + (L_Bp - e_p) * k_Bp) -(r_s + r_p) * ((-L_Ap - e_p) * k_Ap + (L_Bp - e_p) * k_Bp) -kappa_Azp - kappa_Bzp -(r_s + r_p) * (k_Ap + k_Bp) tand(phi_sp)* (r_s + r_p) * (k_Ap + k_Bp) 0; ...
                    0 0 0 0 0 0; ...
                    0 0 0 0 0 0; ...
                    0 0 0 0 0 -k_Azp - k_Bzp];

                % Eq. (28):
                K_c = Upsi_c * sum_sin + R * Upsi_c * R' * sum_cos + Theta_c * sum_sin2 + ...
                    R * Theta_c * R' * sum_cos2 + Xi_c * sum_sc + Psi_c;

                % Eq. (31):
                K_ci = @(i)(Lambda_c(i) * sind(alpha(i)) + R * Lambda_c(i) * cosd(alpha(i)) + Gamma_c(i));

                %% Planet matrices:
                % Eq. (105): relates to planet i. The quantity 2*jdx - 1
                % indicates the sun–planet mesh, and 2*jdx indicates the
                % ring-planet mesh:
                Ki = @(jdx)([k(2 * jdx - 1) * r_p ^ 2 * sind(psi) ^ 2 + k(2 * jdx) * (sind(phi_sp + phi_rp) * ((e_p - c(2 * jdx)) * cosd(psi) - r_p * tand(phi_rp) * sind(psi)) + r_p * sind(psi) * cosd(phi_sp + phi_rp)) ^ 2 + kappa(2 * jdx) * sind(phi_sp + phi_rp) ^ 2 + k_Ap * (-L_Ap - e_p) ^ 2 + k_Bp * (L_Bp - e_p) ^ 2 + kappa_Ap + kappa_Bp k(2 * jdx - 1) * ((e_p - c(2 * jdx - 1)) * cosd(psi) + r_p * tand(phi_sp)* sind(psi)) * r_p * sind(psi) + k(2 * jdx) * (sind(phi_sp + phi_rp) * ((e_p - c(2 * jdx)) * cosd(psi) - r_p * tand(phi_rp) * sind(psi)) + r_p * sind(psi) * cosd(phi_sp + phi_rp)) * (cosd(phi_sp + phi_rp) * ((e_p - c(2 * jdx)) * cosd(psi) - r_p * tand(phi_rp) * sind(psi)) + r_p * sind(psi) * sind(phi_sp + phi_rp)) - kappa(2 * jdx) * cosd(phi_sp + phi_rp) * sind(phi_sp + phi_rp) (k(2 * jdx - 1) * r_p * sind(psi) - k(2 * jdx) * (sind(phi_sp + phi_rp) * ((e_p - c(2 * jdx)) * cosd(psi) - r_p * tand(phi_rp) * sind(psi)) + r_p * sind(psi) * cosd(phi_sp + phi_rp))) * r_p * cosd(psi) -k(2 * jdx - 1) * r_p * sind(psi) * cosd(psi) - k(2 * jdx) * cosd(psi) * cosd(phi_sp + phi_rp) * (sind(phi_sp + phi_rp) * ((e_p - c(2 * jdx)) * cosd(psi) - r_p * tand(phi_rp) * sind(psi)) + r_p * sind(psi) * cosd(phi_sp + phi_rp)) -k(2 * jdx) * cosd(psi) * sind(phi_sp + phi_rp) * (sind(phi_sp + phi_rp) * ((e_p - c(2 * jdx)) * cosd(psi) - r_p * tand(phi_rp) * sind(psi)) + r_p * sind(psi) * cosd(phi_sp + phi_rp)) - (-L_Ap - e_p) * k_Ap - (L_Bp - e_p) * k_Bp k(2 * jdx - 1) * r_p * sind(psi) ^ 2 - k(2 * jdx) * sind(psi) * (sind(phi_sp + phi_rp) * ((e_p - c(2 * jdx)) * cosd(psi) - r_p * tand(phi_rp) * sind(psi)) + r_p * sind(psi) * cosd(phi_sp + phi_rp)); ...
                    0 k(2 * jdx - 1) * ((e_p - c(2 * jdx - 1)) * cosd(psi) + r_p * tand(phi_sp)* sind(psi)) ^ 2 + k(2 * jdx) * (cosd(phi_sp + phi_rp) * ((e_p - c(2 * jdx)) * cosd(psi) - r_p * tand(phi_rp) * sind(psi)) + r_p * sind(psi) * sind(phi_sp + phi_rp)) ^ 2 + kappa(2 * jdx - 1) + kappa(2 * jdx) * cosd(phi_sp + phi_rp) ^ 2 + k_Ap * (-L_Ap - e_p) ^ 2 + k_Bp * (L_Bp - e_p) ^ 2 + kappa_Ap + kappa_Bp (k(2 * jdx - 1) * ((e_p - c(2 * jdx - 1)) * cosd(psi) + r_p * tand(phi_sp)* sind(psi)) - k(2 * jdx) * (cosd(phi_sp + phi_rp) * ((e_p - c(2 * jdx)) * cosd(psi) - r_p * tand(phi_rp) * sind(psi)) + r_p * sind(psi) * sind(phi_sp + phi_rp))) * r_p * cosd(psi) -k(2 * jdx - 1) * ((e_p - c(2 * jdx - 1)) * cosd(psi) + r_p * tand(phi_sp)* sind(psi)) * cosd(psi) - k(2 * jdx) * (cosd(phi_sp + phi_rp) * ((e_p - c(2 * jdx)) * cosd(psi) - r_p * tand(phi_rp) * sind(psi)) + r_p * sind(psi) * sind(phi_sp + phi_rp)) * cosd(psi) * cosd(phi_sp + phi_rp) + (-L_Ap - e_p) * k_Ap + (L_Bp - e_p) * k_Bp -k(2 * jdx) * cosd(psi) * sind(phi_sp + phi_rp) * (cosd(phi_sp + phi_rp) * ((e_p - c(2 * jdx)) * cosd(psi) - r_p * tand(phi_rp) * sind(psi)) + r_p * sind(psi) * sind(phi_sp + phi_rp)) k(2 * jdx - 1) * ((e_p - c(2 * jdx - 1)) * cosd(psi) + r_p * tand(phi_sp)* sind(psi)) * sind(psi) - k(2 * jdx) * (cosd(phi_sp + phi_rp) * ((e_p - c(2 * jdx)) * cosd(psi) - r_p * tand(phi_rp) * sind(psi)) + r_p * sind(psi) * sind(phi_sp + phi_rp)) * sind(psi); ...
                    0 0 (k(2 * jdx - 1) + k(2 * jdx)) * r_p ^ 2 * cosd(psi) ^ 2 + kappa_Azp + kappa_Bzp (k(2 * jdx) * cosd(psi) * cosd(phi_sp + phi_rp) - k(2 * jdx - 1) * cosd(psi)) * r_p * cosd(psi) k(2 * jdx) * cosd(psi) ^ 2 * sind(phi_sp + phi_rp) * r_p (k(2 * jdx - 1) + k(2 * jdx)) * r_p * sind(psi) * cosd(psi); ...
                    0 0 0 k(2 * jdx - 1) * cosd(psi) ^ 2 + k(2 * jdx) * cosd(psi) ^ 2 * cosd(phi_sp + phi_rp) ^ 2 + k_Ap + k_Bp k(2 * jdx) * cosd(psi) ^ 2 * sind(phi_sp + phi_rp) * cosd(phi_sp + phi_rp) -(k(2 * jdx - 1) * cosd(psi) - k(2 * jdx) * cosd(psi) * cosd(phi_sp + phi_rp)) * sind(psi); ...
                    0 0 0 0 cosd(psi) ^ 2 * sind(phi_sp + phi_rp) ^ 2 * k(2 * jdx) + K_Bp + k_Ap k(2 * jdx) * sind(psi) * cosd(psi) * sind(phi_sp + phi_rp); ...
                    0 0 0 0 0 (k(2 * jdx - 1) + k(2 * jdx)) * sind(psi) ^ 2 + k_Azp + k_Bzp]);

                %% Assembling overall matrix:
                np = 6 * stage_idx.N_p;
                K11 = blkdiag(K_s, K_r, K_c);

                K_s_row = zeros(6, np);
                K_r_row = K_s_row;
                K_c_row = K_s_row;
                K_cell = cell(1, stage_idx.N_p);

                for idx = 1:stage_idx.N_p
                    rng = (-5:0) + 6 * idx;
                    K_s_row(:, rng) = K_si(idx);
                    K_r_row(:, rng) = K_ri(idx);
                    K_c_row(:, rng) = K_ci(idx);
                    
                    K_cell{idx} = Ki(idx);
                end

                K12 = [K_s_row; ...
                       K_r_row; ...
                       K_c_row];

                K22 = blkdiag(K_cell{:});

                KK = [K11,  K12; ...
                      K12', K22];

            end

        end
        
        function flag = test(num_planet)
            %TEST
            % +-----------------------------------------------------+-------+--------+--------+--------+--------+---------+
            % | Parameter                                           |  Sun  |  Mesh  | Planet |  Mesh  |  Ring  | Carrier |
            % +-----------------------------------------------------+-------+--------+--------+--------+--------+---------+
            % | Operating pressure angle, phi (deg.)                |       |  21.3  |        |  21.3  |        |         |
            % +-----------------------------------------------------+-------+--------+--------+--------+--------+---------+
            % | Base helix angle, psi (deg.)                        |       |  -28.5 |        |  28.5  |        |         |
            % +-----------------------------------------------------+-------+--------+--------+--------+--------+---------+
            % | Translational mesh stiffness, k (N/m)               |       | 6.19e9 |        | 22.3e9 |        |         |
            % +-----------------------------------------------------+-------+--------+--------+--------+--------+---------+
            % | Tilting mesh stiffness, kappa (N-m/rad)             |       |  643e3 |        | 2.31e6 |        |         |
            % +-----------------------------------------------------+-------+--------+--------+--------+--------+---------+
            % | Center of stiffness, c (mm)                         |       |    0   |        |    0   |        |         |
            % +-----------------------------------------------------+-------+--------+--------+--------+--------+---------+
            % | Base radius, r (mm)                                 |   24  |        |   16   |        |   56   |         |
            % +-----------------------------------------------------+-------+--------+--------+--------+--------+---------+
            % | Center of mass, e (mm)                              |   0   |        |    0   |        |    0   |    0    |
            % +-----------------------------------------------------+-------+--------+--------+--------+--------+---------+
            % | Bearing distance @ A, L^A (mm)                      |  -20  |        |   -20  |        |   -20  |   -20   |
            % +-----------------------------------------------------+-------+--------+--------+--------+--------+---------+
            % | Bearing distance @ B, L^B (mm)                      |   20  |        |   20   |        |   20   |    20   |
            % +-----------------------------------------------------+-------+--------+--------+--------+--------+---------+
            % | Radial bearing stiffness, k^A/B (N/m)               | 500e6 |        |  500e6 |        |  500e6 |  500e6  |
            % +-----------------------------------------------------+-------+--------+--------+--------+--------+---------+
            % | Axial bearing stiffness, K^Az/Bz (N/m)              | 500e6 |        |  500e6 |        |  500e6 |  500e6  |
            % +-----------------------------------------------------+-------+--------+--------+--------+--------+---------+
            % | Tilting bearing stiffness, kappa^A/B (N-m/rad)      |  50e6 |        |   5e6  |        |  50e6  |   50e6  |
            % +-----------------------------------------------------+-------+--------+--------+--------+--------+---------+
            % | Rotational bearing stiffness, kappa^Az/Bz (N-m/rad) |   0   |        |    0   |        |  90e9  |   90e9  |
            % +-----------------------------------------------------+-------+--------+--------+--------+--------+---------+
            % | Mass, m (kg)                                        |  0.3  |        |   0.2  |        | 100e-6 |   0.5   |
            % +-----------------------------------------------------+-------+--------+--------+--------+--------+---------+
            % | Tilting inertia, J^x (kg-m^2)                       |  5e-3 |        |  50e-6 |        |  10e-6 |   4e-3  |
            % +-----------------------------------------------------+-------+--------+--------+--------+--------+---------+
            % | Rotational inertia, J^z (kg-m^2)                    | 10e-3 |        | 100e-6 |        |  20e-6 |   8e-3  |
            % +-----------------------------------------------------+-------+--------+--------+--------+--------+---------+
            %
            % Results:
            % +----+-----------+-----------+
            % |  i | 4 planets | 5 planets |
            % +----+-----------+-----------+
            % |  1 |    953    |    1011   |
            % +----+-----------+-----------+
            % |  2 |    3120   |    3068   |
            % +----+-----------+-----------+
            % |  3 |    3120   |    3068   |
            % +----+-----------+-----------+
            % |  4 |    3251   |    3114   |
            % +----+-----------+-----------+
            % |  5 |    3743   |    3670   |
            % +----+-----------+-----------+
            % |  6 |    5426   |    5184   |
            % +----+-----------+-----------+
            % |  7 |    5426   |    5184   |
            % +----+-----------+-----------+
            % |  8 |    8177   |    8177   |
            % +----+-----------+-----------+
            % |  9 |    8537   |    8177   |
            % +----+-----------+-----------+
            % | 10 |    8537   |    8506   |
            % +----+-----------+-----------+
            %
            
            switch(num_planet)
                case 4
                    fn_ref = [ 953 3120 3120 3251 3743 ...
                              5426 5426 8177 8537 8537];
                    
                case 5
                    fn_ref = [1011 3068 3068 3114 3670 ...
                              5184 5184 8177 8177 8506];
                otherwise
                    error('Range is between 4 and 5.');
            end

            % Defining test cases for:
            % Gear_Set (main):
            import matlab.mock.TestCase;
            gset_TC = TestCase.forInteractiveUse;
            
            % Carrier (aux.):
            carrier_TC = TestCase.forInteractiveUse;
            
            

        end


    end


end
