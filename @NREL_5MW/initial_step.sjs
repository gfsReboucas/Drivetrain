function main()
{
    var model = Spck.currentModel;

    // Apply zero initial condtions:
    var zero_IC = model.findElement("$ST_all_zero");
    zero_IC.copyToModel();

    // Disabling translation and tilting for shafts:
    var switch_AIC = model.findElement('$SVG_loading.$_switch_AIC');
    switch_AIC.str.src = '1';

    // Base excitation: enabled
    var base_flag = model.findElement('$SVG_loading.$_base_excitation_flag');
    base_flag.str.src = '1';
    // 3. Aerodynamic loads: nominal
    var aero_flag = model.findElement('$SVG_loading.$_aero_excitation_flag');
    aero_flag.str.src = '1';
    // 4. Generator loads: nominal
    var gen_flag = model.findElement('$SVG_loading.$_generator_flag');
    gen_flag.str.src = '1';
    // 5. gear force elements: basic
    var mesh_flag = model.findElement('$SVG_loading.$_mesh_flag');
    mesh_flag.str.src = '204';

    model.save();

    // Apply Initial Conditions solver:
    var slv_AIC = Spck.Slv.applyInitialConditions(model);

    // Update state set:
    var state_AIC = Spck.currentModel.findElement("$ST_AIC");
    state_AIC.copyFromModel();

    model.save();

    // Re-enabling translation and tilting for shafts:
    switch_AIC.str.src = '0';
    
    // Calculating static equilibrium state:
    var eq = Spck.Slv.equi(model);

    // Update state set:
    var static_eq = Spck.currentModel.findElement("$ST_static_eq");
    static_eq.copyFromModel();

    // Base excitation: enabled
    base_flag.str.src   = '1';
    // Aerodynamic loads: DLC
    aero_flag.str.src = '2';
    // Generator loads: PI
    gen_flag.str.src = '3';
    // Gear force elements: advanced
    mesh_flag.str.src = '225';

    Spck.currentModel.save();
}