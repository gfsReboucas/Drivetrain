function main()
{
    // Mesh stiffness:
    var RB_mesh_nonlin = new RadioButton;
    RB_mesh_nonlin.text = "Non-linear (FE225)";
    connect(RB_mesh_nonlin, "clicked()",  mesh_flag_nlin);
    connect(RB_mesh_nonlin, "pressed()",  mesh_flag_nlin);
    connect(RB_mesh_nonlin, "released()", mesh_flag_lin);

    var RB_mesh_lin = new RadioButton;
    RB_mesh_lin.text = "Linear (FE204)";
    connect(RB_mesh_lin, "clicked()",  mesh_flag_lin);
    connect(RB_mesh_lin, "pressed()",  mesh_flag_lin);
    connect(RB_mesh_lin, "released()", mesh_flag_nlin);

    var mesh_flag = Spck.currentModel.findElement("$SVG_mesh.$_flag").getDiscrAlternative();

    switch(mesh_flag)
    {
        case 0:
            RB_mesh_nonlin.checked = true;
            RB_mesh_lin.checked = false;
            break;

        case 1:
            RB_mesh_nonlin.checked = false;
            RB_mesh_lin.checked = true;
            break;

        default:
            print("SCRIPT: Error! Value for [mesh_flag] is undefined.");
    }

    var GB_mesh = new GroupBox;
    GB_mesh.title = "Mesh stiffness formulation:";
    GB_mesh.add(RB_mesh_lin);
    GB_mesh.add(RB_mesh_nonlin);
    
    function mesh_flag_lin()
    {
        mesh_flag = 1;
    }

    function mesh_flag_nlin()
    {
        mesh_flag = 0;
    }

    // Generator:
    var RB_gen_speed = new RadioButton;
    RB_gen_speed.text = "Recorded speed";
    connect(RB_gen_speed, "clicked()" , gen_flag_speed);
    connect(RB_gen_speed, "pressed()" , gen_flag_speed);
    connect(RB_gen_speed, "released()", gen_flag_off);

    var RB_gen_PI = new RadioButton;
    RB_gen_PI.text = "PI speed control";
    connect(RB_gen_PI, "clicked()" , gen_flag_PI);
    connect(RB_gen_PI, "pressed()" , gen_flag_PI);
    connect(RB_gen_PI, "released()", gen_flag_off);

    var RB_gen_off = new RadioButton;
    RB_gen_off.text = "Off";
    connect(RB_gen_off, "clicked()" , gen_flag_off);
    connect(RB_gen_off, "pressed()" , gen_flag_off);

    var gen_flag = Spck.currentModel.findElement("$SVG_generator.$_flag").getDiscrAlternative();

    switch(gen_flag)
    {
        case 0:
            RB_gen_speed.checked = true;
            RB_gen_PI.checked = false;
            RB_gen_off.checked = false;
            break;

        case 1:
            RB_gen_speed.checked = false;
            RB_gen_PI.checked = true;
            RB_gen_off.checked = false;
            break;

        case 1:
            RB_gen_speed.checked = false;
            RB_gen_PI.checked = true;
            RB_gen_off.checked = false;
            break;

        case 2:
            RB_gen_speed.checked = false;
            RB_gen_PI.checked = false;
            RB_gen_off.checked = true;
            break;

        default:
            print("SCRIPT: Error! Value for [gen_flag] is undefined.");
    }

    var GB_gen = new GroupBox;
    GB_gen.title = "Generator torque from:";
    GB_gen.add(RB_gen_speed);
    GB_gen.add(RB_gen_PI);
    GB_gen.add(RB_gen_off);

    function gen_flag_speed()
    {
        gen_flag = 0;
    }

    function gen_flag_PI()
    {
        gen_flag = 1;
    }

    function gen_flag_off()
    {
        gen_flag = 2;
    }

    // bed plate:
    var RB_base_on = new RadioButton;
    RB_base_on.text = "On";
    connect(RB_base_on, "clicked()",  base_flag_true);
    connect(RB_base_on, "pressed()",  base_flag_true);
    connect(RB_base_on, "released()", base_flag_false);

    var RB_base_off = new RadioButton;
    RB_base_off.text = "Off";
    connect(RB_base_off, "clicked()",  base_flag_false);
    connect(RB_base_off, "pressed()",  base_flag_false);
    connect(RB_base_off, "released()", base_flag_true);

    var base_flag = Spck.currentModel.findElement("$SVG_base_input.$_flag").getDiscrAlternative();

    switch(base_flag)
    {
        case 0:
            RB_base_on.checked = true;
            RB_base_off.checked = false;
            break;

        case 1:
            RB_base_on.checked = false;
            RB_base_off.checked = true;
            break;

        default:
            print("SCRIPT: Error! Value for [base_flag] is undefined.");
    }

    var GB_base = new GroupBox;
    GB_base.title = "Excitation at bed plate:";
    GB_base.add(RB_base_on);
    GB_base.add(RB_base_off);

    function base_flag_true()
    {
        base_flag = 0;
    }

    function base_flag_false()
    {
        base_flag = 1;
    }

    var diag = new Dialog;
    diag.title = "Set Options for SubVars:";
    diag.add(GB_mesh);
    diag.add(GB_gen);
    diag.add(GB_base);
    diag.visible = true;

    // execute dialog:
    if(diag.exec()) // OK button is pressed
    {
        Spck.currentModel.findElement("$SVG_mesh.$_flag").setDiscrAlternative(mesh_flag);
        Spck.currentModel.findElement("$SVG_generator.$_flag").setDiscrAlternative(gen_flag);
        Spck.currentModel.findElement("$SVG_base_input.$_flag").setDiscrAlternative(base_flag);
        Spck.currentModel.save();

        print("SCRIPT: SubVar options settings were modified.");
    }
    else // Cancel button is pressed
    {
        print("SCRIPT: SubVar options settings were cancelled.");
    }

}