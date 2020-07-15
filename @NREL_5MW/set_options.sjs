function main()
{
    var diag = new Dialog;
    diag.title = "Set Options for subVars";

    // Mesh stiffness:
    var mesh_flag = 0;
    var GB_mesh = new GroupBox;
    GB_mesh.title = "Mesh stiffness formulation:";

    var RB_mesh_lin = new RadioButton;
    RB_mesh_lin.text = "Linear (FE204)";
    RB_mesh_lin.checked = false;
    connect(RB_mesh_lin, "clicked()", mesh_flag_lin);
    connect(RB_mesh_lin, "pressed()", mesh_flag_lin);
    connect(RB_mesh_lin, "released()", mesh_flag_nlin);
    GB_mesh.add(RB_mesh_lin);

    var RB_mesh_nonlin = new RadioButton;
    RB_mesh_nonlin.text = "Non-linear (FE225)";
    RB_mesh_nonlin.checked = true;
    connect(RB_mesh_nonlin, "clicked()", mesh_flag_nlin);
    connect(RB_mesh_nonlin, "pressed()", mesh_flag_nlin);
    GB_mesh.add(RB_mesh_nonlin);
    
    diag.add(GB_mesh);

    function mesh_flag_lin()
    {
        mesh_flag = 1;
    }

    function mesh_flag_nlin()
    {
        mesh_flag = 0;
    }

    // Generator:
    var gen_flag = 1;
     var GB_gen = new GroupBox;
     GB_gen.title = "Generator torque from:";

    var RB_gen_speed = new RadioButton;
    RB_gen_speed.text = "Recorded speed";
    RB_gen_speed.checked = true;
    connect(RB_gen_speed, "clicked()" , gen_flag_speed);
    connect(RB_gen_speed, "pressed()" , gen_flag_speed);
    connect(RB_gen_speed, "released()", gen_flag_off);
    GB_gen.add(RB_gen_speed);
    
    var RB_gen_PI = new RadioButton;
    RB_gen_PI.text = "PI speed control";
    connect(RB_gen_PI, "clicked()" , gen_flag_PI);
    connect(RB_gen_PI, "pressed()" , gen_flag_PI);
    connect(RB_gen_PI, "released()", gen_flag_off);
    GB_gen.add(RB_gen_PI);

    var RB_gen_off = new RadioButton;
    RB_gen_off.text = "Off";
    connect(RB_gen_off, "clicked()" , gen_flag_off);
    connect(RB_gen_off, "pressed()" , gen_flag_off);
    GB_gen.add(RB_gen_off);

    diag.add(GB_gen);

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
    var base_flag = 1;
    var GB_bed = new GroupBox;
    GB_bed.title = "Excitation at bed plate:";

    var RB_bed_on = new RadioButton;
    RB_bed_on.text = "On";
    RB_bed_on.checked = true;
    connect(RB_bed_on, "clicked()" , base_flag_true);
    connect(RB_bed_on, "pressed()" , base_flag_true);
    GB_bed.add(RB_bed_on);

    var RB_bed_off = new RadioButton;
    RB_bed_off.text = "Off";
    connect(RB_bed_off, "clicked()" , base_flag_false);
    connect(RB_bed_off, "pressed()" , base_flag_false);
    GB_bed.add(RB_bed_off);

    diag.add(GB_bed);

    function base_flag_true()
    {
        base_flag = 0;
    }

    function base_flag_false()
    {
        base_flag = 1;
    }

    diag.visible = true;

    // execute dialog:
    if(diag.exec()) // OK button is pressed
    {
        Spck.currentModel.findElement("$SVG_mesh.$_flag").setDiscrAlternative(mesh_flag);
        Spck.currentModel.findElement("$SVG_generator.$_flag").setDiscrAlternative(gen_flag);
        Spck.currentModel.findElement("$SVG_base_input.$_flag").setDiscrAlternative(base_flag);

        print("SCRIPT: SubVar options settings were modified.");
    }
    else // Cancel button is pressed
    {
        print("SCRIPT: SubVar options settings were cancelled.");
    }

}