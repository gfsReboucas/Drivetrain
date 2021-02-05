function main()
{
    // Rack reference profile:
    ref_prof_box = new ComboBox;
    ref_prof_box.label = 'Reference profile, acc. ISO 53:';
    ref_prof_box.itemList = ['A', 'B', 'C', 'D'];

    // normal module:
    var m_n_edit = new NumberEdit;
    m_n_edit.decimals = 3;
    m_n_edit.minimun = 1.0;
    m_n_edit.maximum = 50.0;
    m_n_edit.value = 8.0;
    m_n_edit.label = 'Normal module, m_n in mm:';
    
    var alpha_n_edit = new NumberEdit;
    alpha_n_edit.decimals = 1;
    alpha_n_edit.value = 20.0;
    alpha_n_edit.label = 'Normal pressure angle, alpha_n in deg.:';
    
    var beta_edit = new NumberEdit;
    beta_edit.decimals = 1;
    beta_edit.value = 15.8;
    beta_edit.label = 'Helix angle, beta in deg.:';
    
    var a_w_edit = new NumberEdit;
    a_w_edit.decimals = 1;
    a_w_edit.value = 500.0;
    a_w_edit.label = 'Center distance, a_w in mm:';
    
    var general_box = new GroupBox;
    general_box.title = 'General parameters:';
    general_box.add(ref_prof_box);
    general_box.add(m_n_edit);
    general_box.add(alpha_n_edit);
    general_box.add(beta_edit);
    general_box.add(a_w_edit);

    // pinion:
    var z1_edit = new NumberEdit;
    z1_edit.decimals = 0;
    z1_edit.label = 'Number of teeth, z_1:';
    z1_edit.value = 17;

    var b1_edit = new NumberEdit;
    b1_edit.decimals = 2;
    b1_edit.label = 'Face width, [mm] b_1:';
    b1_edit.value = 100.0;

    var x1_edit = new NumberEdit;
    x1_edit.decimals = 4;
    x1_edit.label = 'Profile shift coeff. x_1:';
    x1_edit.value = 0.145;

    // wheel:
    var z2_edit = new NumberEdit;
    z2_edit.decimals = 0;
    z2_edit.label = 'z_2:';
    z2_edit.value = 103;

    var b2_edit = new NumberEdit;
    b2_edit.decimals = 2;
    b2_edit.label = 'b_2:';
    b2_edit.value = 100.0;

    var x2_edit = new NumberEdit;
    x2_edit.decimals = 4;
    x2_edit.label = 'x_2:';
    x2_edit.value = 0.0;

    var diag = new Dialog;
    diag.title = "Set the gear set's properties:";
    diag.add(general_box);
    
    diag.newTab('Parallel');
    diag.add(z1_edit);
    diag.add(b1_edit);
    diag.add(x1_edit);
    diag.newColumn();
    diag.add(z2_edit);
    diag.add(b2_edit);
    diag.add(x2_edit);

    diag.newTab('Planetary');

    diag.visible = true;

    // execute dialog:
    if(diag.exec()) // OK button is pressed
    {
        Spck.currentModel.save();

        print("SCRIPT: it works.");
    }
    else // Cancel button is pressed
    {
        print("SCRIPT: cancelled operation.");
    }

}