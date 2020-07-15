function main()
{
        var joint = Spck.currentModel.findElement("$J_Bed_Plate"); // get joint

        joint.from = Spck.currentModel.findElement("$M_Isys");          // From Marker
        joint.to   = Spck.currentModel.findElement("$M_Bed_Plate_BRF"); // To Marker

        var flag  = Spck.currentModel.findElement("$SVG_base_input.$_flag").eval(); // get flag
        if(flag == true)
        {
                joint.type.src = 35; // Type
 
                var param = ["X"       , "Y"       , "Z"       , // displacement 
                             "X_Rot"   , "Y_Rot"   , "Z_Rot"   , // rotation
                             "XD"      , "YD"      , "ZD"      , // velocity
                             "X_Rot_D" , "Y_Rot_D" , "Z_Rot_D" , // angular velocity
                             "XDD"     , "YDD"     , "ZDD"     , // acceleration
                             "X_Rot_DD", "Y_Rot_DD", "Z_Rot_DD", // angular acceleration
                             1.0];                               // rotation sequence

                // input functions:
                for(idx = 0; idx < param.length - 1; idx++)
                {
                        param[idx] = "$IFG_bed_plate_motion.$I_" + param[idx];
                        param[idx] = Spck.currentModel.findElement(param[idx]);
                }
        }
        else
        {
                joint.type.src = 0; // Type
                joint.addAttr(".st.vel(1", ") = 0.0"); // Velocity
                
                var param = [0.0, 0.0, 0.0, // rotation about alpha, beta, and gamma [rad]
                             0.0, 0.0, 0.0, // translation in x, y, and z [m]
                             0.0];          // Rotation sequence [-]
        }

        joint.setPar(param);

        joint.addAttr(".2d.pos.x(1,", ") =  320.0");
        joint.addAttr(".2d.pos.y(1,", ") =  -20.0");
        joint.addAttr(".2d.ori  (1,", ") =    0.0");
        joint.addAttr(".2d.paths.from.x(1, 1, 1, 1,", ") =  100.0,  100.0 ,  110.0 ,  120.0, 130.0, 140.0, 283.0,  283.0,  293.0,  303.0");
        joint.addAttr(".2d.paths.from.y(1, 1, 1, 1,", ") =  -10.0,    0.0 ,    0.0 ,    0.0,   0.0,   0.0,   0.0,  -20.0,  -20.0,  -20.0");
        joint.addAttr(".2d.paths.to.x  (1, 1, 1, 1,", ") =  370.0,  354.25,  354.25,  337.0");
        joint.addAttr(".2d.paths.to.y  (1, 1, 1, 1,", ") =  -10.0,  -10.0 ,  -20.0 ,  -20.0");

        Spck.currentModel.findElement("$J_Bed_Plate").src =  joint;
        Spck.currentModel.save();

        print("SCRIPT: joint $J_Bed_Plate updated succesfully.");

}