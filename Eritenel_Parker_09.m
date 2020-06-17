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
        end
    end
    
    
end
