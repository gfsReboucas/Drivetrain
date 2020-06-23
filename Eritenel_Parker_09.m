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
    
    methods
        function obj = Eritenel_Parker_09(DT)
            obj@Dynamic_Formulation(DT);
            
            [obj.M, obj.G] = obj.inertia_matrix();
            obj.K = obj.stiffness_matrix();
            
        end
    end
    
    methods(Static)
        function MM = stage_inertia_matrix(stage_idx)
        end
        
        function KK = stage_stiffness_matrix(stage_idx)
            N1    = stage_idx.N_p - 1;
            
            alpha = (0:N1)*(360.0/stage_idx.N_p);
            sum_sin = sum(sind(alpha));
            sum_cos = sum(cosd(alpha));
            
            R = [0 1 0 0 0 0;
                -1 0 0 0 0 0;
                0 0 1 0 0 0;
                0 0 0 0 1 0;
                0 0 0 -1 0 0;
                0 0 0 0 0 1];
            %%
            Upsi_s = [0 0 k_1 * r_s * ((e_s - c_1) * cosd(psi) - r_s * tand(phi_sp) * sind(psi)) * cosd(psi) 0 0 -k_1 * ((e_s - c_1) * cosd(psi) - r_s * tand(phi_sp) * sind(psi)) * sind(psi);
                0 0 k_1 * r_s ^ 2 * sind(psi) * cosd(psi) 0 0 -k_1 * r_s * sind(psi) ^ 2;
                0 0 0 0 k_1 * r_s * cosd(psi) ^ 2 0;
                0 0 0 0 0 0;
                0 0 0 0 0 0;
                0 0 0 0 0 0];
            
            Theta_s = [kappa_1 + k_1 * ((e_s - c_1) * cosd(psi) - r_s * tand(phi_sp) * sind(psi)) ^ 2 k_1 * ((e_s - c_1) * cosd(psi) - r_s * tand(phi_sp) * sind(psi)) * r_s * sind(psi) 0 0 k_1 * ((e_s - c_1) * cosd(psi) - r_s * tand(phi_sp) * sind(psi)) * cosd(psi) 0;
                0 k_1 * r_s * sind(psi) ^ 2 0 0 k_1 * r_s * sind(psi) * cosd(psi) 0;
                0 0 0 0 0 0;
                0 0 0 0 0 0;
                0 0 0 0 k_1 * cosd(psi) ^ 2 0;
                0 0 0 0 0 0];
            
            Xi_s = [0.2e1 * k_1 * ((e_s - c_1) * cosd(psi) - r_s * tand(phi_sp) * sind(psi)) * r_s * sind(psi) k_1 * (r_s ^ 2 * sind(psi) ^ 2 - ((e_s - c_1) * cosd(psi) - r_s * tand(phi_sp) * sind(psi)) ^ 2) - kappa_1 0 k_1 * ((e_s - c_1) * cosd(psi) - r_s * tand(phi_sp) * sind(psi)) * cosd(psi) k_1 * r_s * sind(psi) * cosd(psi) 0;
                0 -0.2e1 * k_1 * ((e_s - c_1) * cosd(psi) - r_s * tand(phi_sp) * sind(psi)) * r_s * sind(psi) 0 k_1 * r_s * sind(psi) * cosd(psi) -k_1 * ((e_s - c_1) * cosd(psi) - r_s * tand(phi_sp) * sind(psi)) * cosd(psi) 0;
                0 0 0 0 0 0;
                0 0 0 0 k_1 * cosd(psi) ^ 2 0;
                0 0 0 0 0 0;
                0 0 0 0 0 0];
            
            Psi_s = [(-L_As - e_s) ^ 2 * k_As + (L_Bs - e_s) ^ 2 * k_Bs + kappa_As + kappa_BBs 0 0 0 -(-L_As - e_s) * k_As - (L_Bs - e_s) * k_Bs 0;
                0 (-L_As - e_s) ^ 2 * k_As + (L_Bs - e_s) ^ 2 * k_Bs + kappa_As + kappa_BBs 0 (-L_As - e_s) * k_As + (L_Bs - e_s) * k_Bs 0 0;
                0 0 kappa_Azs + kappa_Bzs + k_1 * r_s ^ 2 * cosd(psi) ^ 2 0 0 -k_1 * r_s * sind(psi) * cosd(psi);
                0 0 0 k_As + k_Bs 0 0;
                0 0 0 0 k_As + k_Bs 0;
                0 0 0 0 0 k_Azs + k_Bzs + k_1 * sind(psi) ^ 2];
            
            Lambda_s = [k_1 * ((e_s - c_1) * cosd(psi) - r_s * tand(phi_sp) * sind(psi)) * r_p * sind(psi) k_1 * ((e_s - c_1) * cosd(psi) - r_s * tand(phi_sp) * sind(psi)) * ((e_p - c_1) * cosd(psi) + r_p * tand(phi_sp) * sind(psi)) + kappa_1 k_1 * ((e_s - c_1) * cosd(psi) - r_s * tand(phi_sp) * sind(psi)) * r_p * cosd(psi) -k_1 * ((e_s - c_1) * cosd(psi) - r_s * tand(phi_sp) * sind(psi)) * cosd(psi) 0 k_1 * ((e_s - c_1) * cosd(psi) - r_s * tand(phi_sp) * sind(psi)) * sind(psi);
                k_1 * r_s * r_p * sind(psi) ^ 2 k_1 * r_s * sind(psi) * ((e_p - c_1) * cosd(psi) + r_p * tand(phi_sp) * sind(psi)) k_1 * r_s * r_p * sind(psi) * cosd(psi) -k_1 * r_s * sind(psi) * cosd(psi) 0 k_1 * r_s * sind(psi) ^ 2;
                0 0 0 0 0 0;
                0 0 0 0 0 0;
                k_1 * r_p * sind(psi) * cosd(psi) k_1 * ((e_p - c_1) * cosd(psi) + r_p * tand(phi_sp) * sind(psi)) * cosd(psi) k_1 * r_p * cosd(psi) ^ 2 -k_1 * cosd(psi) ^ 2 0 k_1 * sind(psi) * cosd(psi);
                0 0 0 0 0 0];
            
            %%
            Upsi_r = [0 0 k_2 * r_r * (cosd(phi_sp + phi_rp) * ((e_r - c_2) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * sind(phi_sp + phi_rp)) * cosd(psi) 0 0 k_2 * (cosd(phi_sp + phi_rp) * ((e_r - c_2) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * sind(phi_sp + phi_rp)) * sind(psi);
                0 0 k_2 * r_r * (sind(phi_sp + phi_rp) * ((e_r - c_2) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * cosd(phi_sp + phi_rp)) * cosd(psi) 0 0 k_2 * (sind(phi_sp + phi_rp) * ((e_r - c_2) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * cosd(phi_sp + phi_rp)) * sind(psi);
                0 0 0 -k_2 * r_r * cosd(psi) ^ 2 * sind(phi_sp + phi_rp) k_2 * r_r * cosd(psi) ^ 2 * cosd(phi_sp + phi_rp) 0;
                0 0 0 0 0 -k_2 * sind(psi) * cosd(psi) * sind(phi_sp + phi_rp);
                0 0 0 0 0 k_2 * cosd(psi) * cosd(phi_sp + phi_rp) * sind(psi);
                0 0 0 0 0 0;];
            
            Theta_r = [kappa_2 * cosd(phi_sp + phi_rp) ^ 2 + k_2 * (cosd(phi_sp + phi_rp) * ((e_r - c_2) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * sind(phi_sp + phi_rp)) ^ 2 kappa_2 * cosd(phi_sp + phi_rp) * sind(phi_sp + phi_rp) + k_2 * (cosd(phi_sp + phi_rp) * ((e_r - c_2) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * sind(phi_sp + phi_rp)) * (sind(phi_sp + phi_rp) * ((e_r - c_2) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * cosd(phi_sp + phi_rp)) 0 -k_2 * (cosd(phi_sp + phi_rp) * ((e_r - c_2) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * sind(phi_sp + phi_rp)) * cosd(psi) * sind(phi_sp + phi_rp) k_2 * (cosd(phi_sp + phi_rp) * ((e_r - c_2) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * sind(phi_sp + phi_rp)) * cosd(psi) * cosd(phi_sp + phi_rp) 0;
                0 kappa_2 * sind(phi_sp + phi_rp) ^ 2 + k_2 * (sind(phi_sp + phi_rp) * ((e_r - c_2) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * cosd(phi_sp + phi_rp)) ^ 2 0 -k_2 * (sind(phi_sp + phi_rp) * ((e_r - c_2) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * cosd(phi_sp + phi_rp)) * cosd(psi) * sind(phi_sp + phi_rp) k_2 * (sind(phi_sp + phi_rp) * ((e_r - c_2) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * cosd(phi_sp + phi_rp)) * cosd(psi) * cosd(phi_sp + phi_rp) 0;
                0 0 0 0 0 0;
                0 0 0 k_2 * cosd(psi) ^ 2 * sind(phi_sp + phi_rp) ^ 2 -k_2 * cosd(psi) ^ 2 * sind(phi_sp + phi_rp) * cosd(phi_sp + phi_rp) 0;
                0 0 0 0 k_2 * cosd(psi) ^ 2 * cosd(phi_sp + phi_rp) ^ 2 0;
                0 0 0 0 0 0;];
            
            Xi_r = [0.2e1 * kappa_2 * cosd(phi_sp + phi_rp) * sind(phi_sp + phi_rp) + 0.2e1 * k_2 * (cosd(phi_sp + phi_rp) * ((e_r - c_2) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * sind(phi_sp + phi_rp)) * (sind(phi_sp + phi_rp) * ((e_r - c_2) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * cosd(phi_sp + phi_rp)) k_2 * ((sind(phi_sp + phi_rp) * ((e_r - c_2) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * cosd(phi_sp + phi_rp)) ^ 2 - (cosd(phi_sp + phi_rp) * ((e_r - c_2) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * sind(phi_sp + phi_rp)) ^ 2) + kappa_2 * (sind(phi_sp + phi_rp) ^ 2 - cosd(phi_sp + phi_rp) ^ 2) 0 k_2 * ((cosd(phi_sp + phi_rp) * ((e_r - c_2) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * sind(phi_sp + phi_rp)) * cosd(psi) * cosd(phi_sp + phi_rp) - (sind(phi_sp + phi_rp) * ((e_r - c_2) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * cosd(phi_sp + phi_rp)) * cosd(psi) * sind(phi_sp + phi_rp)) k_2 * ((sind(phi_sp + phi_rp) * ((e_r - c_2) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * cosd(phi_sp + phi_rp)) * cosd(psi) * cosd(phi_sp + phi_rp) - (cosd(phi_sp + phi_rp) * ((e_r - c_2) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * sind(phi_sp + phi_rp)) * cosd(psi) * sind(phi_sp + phi_rp)) 0;
                0 -0.2e1 * kappa_2 * cosd(phi_sp + phi_rp) * sind(phi_sp + phi_rp) - 0.2e1 * k_2 * (cosd(phi_sp + phi_rp) * ((e_r - c_2) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * sind(phi_sp + phi_rp)) * (sind(phi_sp + phi_rp) * ((e_r - c_2) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * cosd(phi_sp + phi_rp)) 0 k_2 * ((sind(phi_sp + phi_rp) * ((e_r - c_2) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * cosd(phi_sp + phi_rp)) * cosd(psi) * cosd(phi_sp + phi_rp) - (cosd(phi_sp + phi_rp) * ((e_r - c_2) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * sind(phi_sp + phi_rp)) * cosd(psi) * sind(phi_sp + phi_rp)) -k_2 * ((cosd(phi_sp + phi_rp) * ((e_r - c_2) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * sind(phi_sp + phi_rp)) * cosd(psi) * cosd(phi_sp + phi_rp) - (sind(phi_sp + phi_rp) * ((e_r - c_2) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * cosd(phi_sp + phi_rp)) * cosd(psi) * sind(phi_sp + phi_rp)) 0;
                0 0 0 0 0 0;
                0 0 0 -0.2e1 * k_2 * cosd(psi) ^ 2 * sind(phi_sp + phi_rp) * cosd(phi_sp + phi_rp) k_2 * (-cosd(psi) ^ 2 * sind(phi_sp + phi_rp) ^ 2 + cosd(psi) ^ 2 * cosd(phi_sp + phi_rp) ^ 2) 0;
                0 0 0 0 0.2e1 * k_2 * cosd(psi) ^ 2 * sind(phi_sp + phi_rp) * cosd(phi_sp + phi_rp) 0;
                0 0 0 0 0 0];
            
            Psi_r = [k_Ar * (-L_Ar - e_r) ^ 2 + k_Br * (L_Br - e_r) ^ 2 + kappa_Ar + kappa_Br 0   0   0   -(-L_Ar - e_r) * k_Ar - (L_Br - e_r) * k_Br   0;
                0  k_Ar * (-L_Ar - e_r) ^ 2 + k_Br * (L_Br - e_r) ^ 2 + kappa_Ar + kappa_Br  0   (-L_Ar - e_r) * k_Ar + (L_Br - e_r) * k_Br    0   0;
                0  0   kappa_Azr + kappa_Bzr + k_2 * r_r ^ 2 * cosd(psi) ^ 2    0   0   k_2 * r_r * cosd(psi) * sind(psi);
                0  (-L_Ar - e_r) * k_Ar + (L_Br - e_r) * k_Br    0   k_Ar + k_Br 0   0;
                -(-L_Ar - e_r) * k_Ar - (L_Br - e_r) * k_Br  0   0   0   k_Ar + k_Br 0;
                0  0   k_2 * r_r * cosd(psi) * sind(psi)   0   0   k_Azr + k_Bzr + k_2 * sind(psi) ^ 2];
            
            Lambda_r = [k_2 * (cosd(phi_sp + phi_rp) * ((e_r - c_2) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * sind(phi_sp + phi_rp)) * (sind(phi_sp + phi_rp) * ((e_p - c_2) * cosd(psi) - r_p * tand(phi_rp) * sind(psi)) + r_p * sind(psi) * cosd(phi_sp + phi_rp)) - kappa_2 * cosd(phi_sp + phi_rp) * sind(phi_sp + phi_rp) k_2 * (cosd(phi_sp + phi_rp) * ((e_r - c_2) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * sind(phi_sp + phi_rp)) * (cosd(phi_sp + phi_rp) * ((e_p - c_2) * cosd(psi) - r_p * tand(phi_rp) * sind(psi)) + r_p * sind(psi) * sind(phi_sp + phi_rp)) - kappa_2 * cosd(phi_sp + phi_rp) ^ 2 -k_2 * r_p * (cosd(phi_sp + phi_rp) * ((e_r - c_2) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * sind(phi_sp + phi_rp)) * cosd(psi) -k_2 * (cosd(phi_sp + phi_rp) * ((e_r - c_2) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * sind(phi_sp + phi_rp)) * cosd(psi) * cosd(phi_sp + phi_rp) -k_2 * (cosd(phi_sp + phi_rp) * ((e_r - c_2) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * sind(phi_sp + phi_rp)) * cosd(psi) * sind(phi_sp + phi_rp) -k_2 * (cosd(phi_sp + phi_rp) * ((e_r - c_2) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * sind(phi_sp + phi_rp)) * sind(psi);
                k_2 * (sind(phi_sp + phi_rp) * ((e_r - c_2) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * cosd(phi_sp + phi_rp)) * (sind(phi_sp + phi_rp) * ((e_p - c_2) * cosd(psi) - r_p * tand(phi_rp) * sind(psi)) + r_p * sind(psi) * cosd(phi_sp + phi_rp)) - kappa_2 * sind(phi_sp + phi_rp) ^ 2 k_2 * (sind(phi_sp + phi_rp) * ((e_r - c_2) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * cosd(phi_sp + phi_rp)) * (cosd(phi_sp + phi_rp) * ((e_p - c_2) * cosd(psi) - r_p * tand(phi_rp) * sind(psi)) + r_p * sind(psi) * sind(phi_sp + phi_rp)) + kappa_2 * cosd(phi_sp + phi_rp) * sind(phi_sp + phi_rp) -k_2 * (sind(phi_sp + phi_rp) * ((e_r - c_2) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * cosd(phi_sp + phi_rp)) * r_p * cosd(psi) -k_2 * (sind(phi_sp + phi_rp) * ((e_r - c_2) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * cosd(phi_sp + phi_rp)) * cosd(psi) * cosd(phi_sp + phi_rp) -k_2 * (sind(phi_sp + phi_rp) * ((e_r - c_2) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * cosd(phi_sp + phi_rp)) * cosd(psi) * sind(phi_sp + phi_rp) -k_2 * (sind(phi_sp + phi_rp) * ((e_r - c_2) * cosd(psi) - r_r * tand(phi_rp) * sind(psi)) - r_r * sind(psi) * cosd(phi_sp + phi_rp)) * sind(psi);
                0 0 0 0 0 0;
                -k_2 * cosd(psi) * sind(phi_sp + phi_rp) * (sind(phi_sp + phi_rp) * ((e_p - c_2) * cosd(psi) - r_p * tand(phi_rp) * sind(psi)) + r_p * sind(psi) * cosd(phi_sp + phi_rp)) -k_2 * cosd(psi) * sind(phi_sp + phi_rp) * (cosd(phi_sp + phi_rp) * ((e_p - c_2) * cosd(psi) - r_p * tand(phi_rp) * sind(psi)) + r_p * sind(psi) * sind(phi_sp + phi_rp)) -k_2 * cosd(psi) ^ 2 * sind(phi_sp + phi_rp) * r_p k_2 * cosd(psi) ^ 2 * sind(phi_sp + phi_rp) * cosd(phi_sp + phi_rp) k_2 * cosd(psi) ^ 2 * sind(phi_sp + phi_rp) ^ 2 k_2 * sind(psi) * cosd(psi) * sind(phi_sp + phi_rp);
                k_2 * cosd(psi) * cosd(phi_sp + phi_rp) * (sind(phi_sp + phi_rp) * ((e_p - c_2) * cosd(psi) - r_p * tand(phi_rp) * sind(psi)) + r_p * sind(psi) * cosd(phi_sp + phi_rp)) k_2 * (cosd(phi_sp + phi_rp) * ((e_p - c_2) * cosd(psi) - r_p * tand(phi_rp) * sind(psi)) + r_p * sind(psi) * sind(phi_sp + phi_rp)) * cosd(psi) * cosd(phi_sp + phi_rp) k_2 * cosd(psi) ^ 2 * cosd(phi_sp + phi_rp) * r_p -k_2 * cosd(psi) ^ 2 * cosd(phi_sp + phi_rp) ^ 2 -k_2 * cosd(psi) ^ 2 * sind(phi_sp + phi_rp) * cosd(phi_sp + phi_rp) -k_2 * cosd(psi) * cosd(phi_sp + phi_rp) * sind(psi);
                0 0 0 0 0 0];
            
            %%
            Upsi_c = [0 0 -(r_s + r_p) * ((-L_Ap - e_s) * k_Ap + (L_Bp - e_s) * k_Bp) 0 0 tand(phi_sp) * (r_s + r_p) * (k_Azp + k_Bzp);
                0 0 -tand(phi_sp) * (r_s + r_p) * ((-L_Ap - e_s) * k_Ap + (L_Bp - e_s) * k_Bp) 0 0 -(r_s + r_p) * (k_Azp + k_Bzp);
                0 0 0 -tand(phi_sp) * (r_s + r_p) * (k_Ap + k_Bp) (r_s + r_p) * (k_Ap + k_Bp) 0;
                0 0 0 0 0 0;
                0 0 0 0 0 0;
                0 0 0 0 0 0];
            
            Theta_c = [k_Ap * (-L_Ap - e_s) ^ 2 + k_Bp * (L_Bp - e_s) ^ 2 + (k_Azp + k_Bzp) * tand(phi_sp) ^ 2 * (r_s + r_p) ^ 2 + kappa_Ap + kappa_Bp -tand(phi_sp) * (r_s + r_p) ^ 2 * (k_Azp + k_Bzp) 0 0 -(-L_Ap - e_s) * k_Ap - (L_Bp - e_s) * k_Bp 0;
                0 k_Ap * (-L_Ap - e_s) ^ 2 + k_Bp * (L_Bp - e_s) ^ 2 + (k_Azp + k_Bzp) * (r_s + r_p) ^ 2 + kappa_Ap + kappa_Bp 0 (-L_Ap - e_s) * k_Ap + (L_Bp - e_s) * k_Bp 0 0;
                0 0 0 0 0 0;
                0 0 0 k_Ap + k_Bp 0 0;
                0 0 0 0 k_Ap + k_Bp 0;
                0 0 0 0 0 0];
            
            Xi_c = [-0.2e1 * tand(phi_sp) * (r_s + r_p) ^ 2 * (k_Azp + k_Bzp) (-tand(phi_sp) ^ 2 * (r_s + r_p) ^ 2 + (r_s + r_p) ^ 2) * (k_Azp + k_Bzp) 0 0 0 0;
                0 0.2e1 * tand(phi_sp) * (r_s + r_p) ^ 2 * (k_Azp + k_Bzp) 0 0 0 0;
                0 0 0 0 0 0;
                0 0 0 0 0 0;
                0 0 0 0 0 0;
                0 0 0 0 0 0];
            
            Psi_c = [(-L_Ac - e_c) ^ 2 * k_Ac + (L_Bc - e_c) ^ 2 * k_Bc + kappa_Ac + kappa_Bc 0 0 0 -(-L_Ac - e_c) * k_Ac - (L_Bc - e_c) * k_Bc 0;
                0 (-L_Ac - e_c) ^ 2 * k_Ac + (L_Bc - e_c) ^ 2 * k_Bc + kappa_Ac + kappa_Bc 0 (-L_Ac - e_c) * k_Ac + (L_Bc - e_c) * k_Bc 0 0;
                0 0 kappa(Az,c) + kappa(Bz,c) + kappa_Azp + kappa_Bzp + (tand(phi_sp) ^ 2 * (r_s + r_p) ^ 2 + (r_s + r_p) ^ 2) * (k_Ap + k_Bp) 0 0 0;
                0 0 0 k_Ac + k_Bc 0 0;
                0 0 0 0 k_Ac + k_Bc 0;
                0 0 0 0 0 k_Azc + k_Bzc + k_Azp + k_Bzp];
            
            Lambda_c = [0 (-L_Ap - e_s) * (-L_Ap - e_p) * k_Ap + (L_Bp - e_s) * (L_Bp - e_p) * k_Bp + kappa_Ap + kappa_Bp 0 (-L_Ap - e_s) * k_Ap + (L_Bp - e_s) * k_Bp 0 -tand(phi_sp) * (r_s + r_p) * (k_Azp + k_Bzp);
                -(-L_Ap - e_s) * (-L_Ap - e_p) * k_Ap - (L_Bp - e_s) * (L_Bp - e_p) * k_Bp - kappa_Ap - kappa_Bp 0 0 0 (-L_Ap - e_s) * k_Ap + (L_Bp - e_s) * k_Bp (r_s + r_p) * (k_Azp + k_Bzp);
                tand(phi_sp) * (r_s + r_p) * ((-L_Ap - e_p) * k_Ap + (L_Bp - e_p) * k_Bp) -(r_s + r_p) * ((-L_Ap - e_p) * k_Ap + (L_Bp - e_p) * k_Bp) -kappa_Azp - kappa_Bzp -(r_s + r_p) * (k_Ap + k_Bp) -tand(phi_sp) * (r_s + r_p) * (k_Ap + k_Bp) 0;
                -(-L_Ap - e_p) * k_Ap - (L_Bp - e_p) * k_Bp 0 0 0 k_Ap + k_Bp 0;
                0 -(-L_Ap - e_p) * k_Ap - (L_Bp - e_p) * k_Bp 0 -k_Ap - k_Bp 0 0;
                0 0 0 0 0 -k_Azp - k_Bzp];

            Ki = [k_1 * r_p ^ 2 * sind(psi) ^ 2 + k_2 * (sind(phi_sp + phi_rp) * ((e_p - c_2) * cosd(psi) - r_p * tand(phi_rp) * sind(psi)) + r_p * sind(psi) * cosd(phi_sp + phi_rp)) ^ 2 + kappa_2 * sind(phi_sp + phi_rp) ^ 2 + k_Ap * (-L_Ap - e_p) ^ 2 + k_Bp * (L_Bp - e_p) ^ 2 + kappa_Ap + kappa_Bp k_1 * ((e_p - c_1) * cosd(psi) + r_p * tand(phi_sp) * sind(psi)) * r_p * sind(psi) + k_2 * (sind(phi_sp + phi_rp) * ((e_p - c_2) * cosd(psi) - r_p * tand(phi_rp) * sind(psi)) + r_p * sind(psi) * cosd(phi_sp + phi_rp)) * (cosd(phi_sp + phi_rp) * ((e_p - c_2) * cosd(psi) - r_p * tand(phi_rp) * sind(psi)) + r_p * sind(psi) * sind(phi_sp + phi_rp)) - kappa_2 * cosd(phi_sp + phi_rp) * sind(phi_sp + phi_rp) (k_1 * r_p * sind(psi) - k_2 * (sind(phi_sp + phi_rp) * ((e_p - c_2) * cosd(psi) - r_p * tand(phi_rp) * sind(psi)) + r_p * sind(psi) * cosd(phi_sp + phi_rp))) * r_p * cosd(psi) -k_1 * r_p * sind(psi) * cosd(psi) - k_2 * cosd(psi) * cosd(phi_sp + phi_rp) * (sind(phi_sp + phi_rp) * ((e_p - c_2) * cosd(psi) - r_p * tand(phi_rp) * sind(psi)) + r_p * sind(psi) * cosd(phi_sp + phi_rp)) -k_2 * cosd(psi) * sind(phi_sp + phi_rp) * (sind(phi_sp + phi_rp) * ((e_p - c_2) * cosd(psi) - r_p * tand(phi_rp) * sind(psi)) + r_p * sind(psi) * cosd(phi_sp + phi_rp)) - (-L_Ap - e_p) * k_Ap - (L_Bp - e_p) * k_Bp k_1 * r_p * sind(psi) ^ 2 - k_2 * sind(psi) * (sind(phi_sp + phi_rp) * ((e_p - c_2) * cosd(psi) - r_p * tand(phi_rp) * sind(psi)) + r_p * sind(psi) * cosd(phi_sp + phi_rp));
                  0 k_1 * ((e_p - c_1) * cosd(psi) + r_p * tand(phi_sp) * sind(psi)) ^ 2 + k_2 * (cosd(phi_sp + phi_rp) * ((e_p - c_2) * cosd(psi) - r_p * tand(phi_rp) * sind(psi)) + r_p * sind(psi) * sind(phi_sp + phi_rp)) ^ 2 + kappa_1 + kappa_2 * cosd(phi_sp + phi_rp) ^ 2 + k_Ap * (-L_Ap - e_p) ^ 2 + k_Bp * (L_Bp - e_p) ^ 2 + kappa_Ap + kappa_Bp (k_1 * ((e_p - c_1) * cosd(psi) + r_p * tand(phi_sp) * sind(psi)) - k_2 * (cosd(phi_sp + phi_rp) * ((e_p - c_2) * cosd(psi) - r_p * tand(phi_rp) * sind(psi)) + r_p * sind(psi) * sind(phi_sp + phi_rp))) * r_p * cosd(psi) -k_1 * ((e_p - c_1) * cosd(psi) + r_p * tand(phi_sp) * sind(psi)) * cosd(psi) - k_2 * (cosd(phi_sp + phi_rp) * ((e_p - c_2) * cosd(psi) - r_p * tand(phi_rp) * sind(psi)) + r_p * sind(psi) * sind(phi_sp + phi_rp)) * cosd(psi) * cosd(phi_sp + phi_rp) + (-L_Ap - e_p) * k_Ap + (L_Bp - e_p) * k_Bp -k_2 * cosd(psi) * sind(phi_sp + phi_rp) * (cosd(phi_sp + phi_rp) * ((e_p - c_2) * cosd(psi) - r_p * tand(phi_rp) * sind(psi)) + r_p * sind(psi) * sind(phi_sp + phi_rp)) k_1 * ((e_p - c_1) * cosd(psi) + r_p * tand(phi_sp) * sind(psi)) * sind(psi) - k_2 * (cosd(phi_sp + phi_rp) * ((e_p - c_2) * cosd(psi) - r_p * tand(phi_rp) * sind(psi)) + r_p * sind(psi) * sind(phi_sp + phi_rp)) * sind(psi);
                  0 0 (k_1 + k_2) * r_p ^ 2 * cosd(psi) ^ 2 + kappa(Az,p) + kappa(Bz,p) (k_2 * cosd(psi) * cosd(phi_sp + phi_rp) - k_1 * cosd(psi)) * r_p * cosd(psi) k_2 * cosd(psi) ^ 2 * sind(phi_sp + phi_rp) * r_p (k_1 + k_2) * r_p * sind(psi) * cosd(psi);
                  0 0 0 k_1 * cosd(psi) ^ 2 + k_2 * cosd(psi) ^ 2 * cosd(phi_sp + phi_rp) ^ 2 + k_Ap + k_Bp k_2 * cosd(psi) ^ 2 * sind(phi_sp + phi_rp) * cosd(phi_sp + phi_rp) -(k_1 * cosd(psi) - k_2 * cosd(psi) * cosd(phi_sp + phi_rp)) * sind(psi);
                  0 0 0 0 cosd(psi) ^ 2 * sind(phi_sp + phi_rp) ^ 2 * k_2 + K_Bp + k_Ap k_2 * sind(psi) * cosd(psi) * sind(phi_sp + phi_rp);
                  0 0 0 0 0 (k_1 + k_2) * sind(psi) ^ 2 + k_Azp + k_Bzp];
            
        end
        
        
    end
    
    
end
