function main()
{
    print('SCRIPT:\tSTART initial_step.sjs')
    var model = Spck.currentModel;

    // Disabling translation and tilting for shafts:
    print('SCRIPT:\t1. Disabling shaft translation and tilting...')
    var switch_AIC = model.findElement('$SVG_loading.$_switch_AIC');
    switch_AIC.str.src = '1';

    // Apply zero initial condtions:
    print('SCRIPT:\t2. Setting states to ZERO...');
    var zero_IC = model.findElement("$ST_all_zero");
    zero_IC.copyToModel();

    // Gear force elements: basic
    var mesh_force = 204;
    print('SCRIPT:\t2.1. Setting gear mesh force to type ' + mesh_force + '...');
    var mesh_flag = model.findElement('$SVG_loading.$_mesh_flag');
    mesh_flag.str.src = '' + mesh_force + '';

    print('SCRIPT:\t2.2. Saving model...');
    model.save();

    // Apply Initial Conditions solver:
    print('SCRIPT:\t3. Running Assemble System solver...');
    Spck.Slv.assmbl(model);
//    Spck.Slv.applyInitialConditions(model);

    // Update state set:
    print('SCRIPT:\t3.1. Saving model with updated state sets from solver...');
    var state_AIC = model.findElement("$ST_AIC");
    state_AIC.copyFromModel();

    model.save();

    // Re-enabling translation and tilting for shafts:
    print('SCRIPT:\t4. Re-enabling shaft translation and tilting...')
    switch_AIC.str.src = '0';

    // Calculating static equilibrium:
    print('SCRIPT:\t5. Calculating static equilibrium using:');
    print('SCRIPT:\t5.1. Time integration method...');
    model.slv.active.src = "$SLV_static_eq_time_integ";
    model.save();

    Spck.Slv.equi(model);

    // Update state set:
    print('SCRIPT:\t5.1.1. Updating state sets with result from solver (static equilibrium via time integration)...');
    var static_eq = model.findElement("$ST_static_eq");
    
    var output_path = model.origin;
    var pos = output_path.lastIndexOf('/');
    output_path = output_path.slice(0, pos);
    output_path = output_path + '/output/' + model.name;

    static_eq.importStateFile(output_path + '.spckst');
//    static_eq.importStateFile('M:/Documents/coding/Drivetrain_2/@NREL_5MW/output/' + model.name + '.spckst');
    static_eq.copyToModel();

//    print('SCRIPT:\t5.1.2. Disabling shaft translation and tilting...')
//    switch_AIC.str.src = '1';

    print('SCRIPT:\t5.2. Newton method...');
    model.slv.active.src = "$SLV_static_eq_Newton_acc";
    model.save();
    Spck.Slv.equi(model);

//    print('SCRIPT:\t5.2.1. Re-enabling shaft translation and tilting...')
//    switch_AIC.str.src = '0';

    // Update state set:
    print('SCRIPT:\t5.2.2. Updating state sets with results from solver (static equilibrium via Newton method and convergence w.r.t. acceleration)...');
    static_eq.importStateFile(output_path + '.spckst');
    static_eq.copyToModel();

    // Gear force elements:
//    mesh_force = 225;
//    print('SCRIPT:\t5.3. Setting gear mesh force to type ' + mesh_force + '...');
//    mesh_flag.str.src = '' + mesh_force + '';

print('SCRIPT:\t6. Saving model...');
    model.slv.active.src = "$SLV_SolverSettings";
    model.save();
    Spck.currentModel.save();

    print('SCRIPT:\tEND initial_step.sjs')
}
