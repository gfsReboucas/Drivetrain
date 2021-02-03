function main()
{
    print('SCRIPT:\tSTART initial_step.sjs')
    var model = Spck.currentModel;

    // Disabling translation and tilting for shafts:
    print('SCRIPT:\t1. Disabling shaft translation and tilting...')
    var switch_AIC = Spck.currentModel.findElement('$SVG_loading.$_switch_AIC');
    switch_AIC.str.src = '1';

    // Apply zero initial condtions:
    print('SCRIPT:\t2. Setting states to ZERO...');
    var zero_IC = Spck.currentModel.findElement("$ST_all_zero");
    zero_IC.copyToModel();

    // Gear force elements: basic
    print('SCRIPT:\t2.1. Setting gear mesh force to type 204 (basic)...');
    var mesh_flag = Spck.currentModel.findElement('$SVG_loading.$_mesh_flag');
    mesh_flag.str.src = '204';

    print('SCRIPT:\t2.2. Saving model...');
    Spck.currentModel.save();

    // Apply Initial Conditions solver:
    print('SCRIPT:\t3. Running Assemble System solver...');
    Spck.Slv.assmbl(Spck.currentModel);

    // Update state set:
    print('SCRIPT:\t3.1. Saving model with updated state sets from solver...');
    var state_AIC = Spck.currentModel.findElement("$ST_AIC");
    state_AIC.copyFromModel();

    Spck.currentModel.save();

    // Re-enabling translation and tilting for shafts:
    print('SCRIPT:\t4. Re-enabling shaft translation and tilting...');
    switch_AIC.str.src = '0';

    // Calculating static equilibrium:
    print('SCRIPT:\t5. Calculating static equilibrium using:');
    print('SCRIPT:\t5.1. Time integration method...');
    Spck.currentModel.slv.active.src = "$SLV_static_eq_time_integ";
    Spck.currentModel.save();

    Spck.Slv.equi(Spck.currentModel);

    // Update state set:
    print('SCRIPT:\t5.1.1. Updating state sets with result from solver (static equilibrium via time integration)...');
    
    var output_path = Spck.currentModel.origin;
    var pos = output_path.lastIndexOf('/');
    output_path = output_path.slice(0, pos);
    output_path = output_path + '/output/' + Spck.currentModel.name;

    var static_eq = Spck.currentModel.findElement("$ST_static_eq");
    static_eq.importStateFile(output_path + '.spckst');
//    static_eq.importStateFile('M:/Documents/coding/Drivetrain_2/@NREL_5MW/output/' + Spck.currentModel.name + '.spckst');
    static_eq.copyToModel();

//    print('SCRIPT:\t5.1.2. Disabling shaft translation and tilting...')
//    switch_AIC.str.src = '1';

    print('SCRIPT:\t5.2. Newton method...');
    Spck.currentModel.slv.active.src = "$SLV_static_eq_Newton_acc";
    Spck.currentModel.save();
    Spck.Slv.equi(Spck.currentModel);

//    print('SCRIPT:\t5.2.1. Re-enabling shaft translation and tilting...')
//    switch_AIC.str.src = '0';

    // Update state set:
    print('SCRIPT:\t5.2.2. Updating state sets with results from solver (static equilibrium via Newton method and convergence w.r.t. acceleration)...');
    static_eq.importStateFile(output_path + '.spckst');
    static_eq.copyToModel();

    print('SCRIPT:\t5.3. Setting gear mesh force to type 225 (advanced)...');
    mesh_flag.str.src = '225';
    
    print('SCRIPT:\t6. Saving model...');
    Spck.currentModel.slv.active.src = "$SLV_SolverSettings";
    Spck.currentModel.save();

    print('SCRIPT:\tEND initial_step.sjs')
}
