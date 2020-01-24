using System;
using System.Linq;

namespace KISSsoft_CSharp
{
    class Gear_Set : Gear
    {
        // attributes
        protected string configuration;
        protected double a_w;
        protected int N_p;
        protected KISSsoftCOM.CKISSsoft KS;

        public Gear_Set(string config, int np, double aw, double mn, double alp, string tp, int[] zz, double bb, double[] xx, double bet) : base(mn, alp, tp, zz, bb, xx, bet)
        {
            this.configuration = config;
            this.a_w = aw;

            KISSsoftCOM.CKISSsoft ks = new KISSsoftCOM.CKISSsoft();

            int n_elem = 0;

            if (this.configuration == "parallel")
            {
                ks.GetModule("Z012", false);

                n_elem = 2; // number of elements
                this.N_p = 1; // number of planets
            }
            else if(this.configuration == "planetary")
            {
                ks.GetModule("Z014", false);

                n_elem = 3; // number of elements
                this.N_p = np; // number of planets
            }

            ks.SetVar("ZS.Geo.mn"  , this.m_n.ToString("G5"));
            ks.SetVar("ZP[0].a"    , this.a_w.ToString("G5"));
            ks.SetVar("ZS.Geo.alfn", this.alpha_n.ToString("G5"));
            ks.SetVar("ZS.Geo.beta", this.beta.ToString("G5"));

            for (int idx = 0; idx < n_elem; idx++)
            {
                ks.SetVar($"ZR[{idx}].z", this.z[idx].ToString("G5"));
                ks.SetVar($"ZR[{idx}].x.nul", this.x[idx].ToString("G5"));
                ks.SetVar($"ZR[{idx}].b", this.b.ToString("G5"));
                ks.SetVar($"ZR[{idx}].Vqual", "6"); // Accuracy grade Q
                ks.SetVar($"ZR[{idx}].mat.bez", "18CrNiMo7-6"); // gear material
                ks.SetVar($"ZR[{idx}].RAH", "0.8"); // Roughness average value DS, flank (μm)
            }

            ks.SetVar("ZR[0].Tool.RefProfile.name", "1.25 / 0.38 / 1.0 ISO 53:1998 Profil A");
            ks.SetVar("ZS.Oil.SchmierTyp", "Oil: ISO-VG 220");
            ks.SetVar("RechSt.ISO6336_2017", "0");
            ks.SetVar("ZR[0].mat.comment", "ISO 6336-5 Figure 9/10 (MQ), Core hardness >=25HRC Jominy J=12mm<HRC28");

            this.KS = ks;
        }

        public double[] structural_integrity(double P_r, double n_r, double K_A, double H)
        {
            this.KS.SetVar("ZS.Pnominal", P_r.ToString("G5"));
            this.KS.SetVar("ZS.nnominal", n_r.ToString("G5"));
            this.KS.SetVar("ZS.KA", K_A.ToString("G5"));
            this.KS.SetVar("ZS.H", H.ToString("G5"));

            this.KS.Calculate();

            int n_elem = 0;

            if (this.configuration == "parallel")
                n_elem = 2;
            else if (this.configuration == "planetary")
                n_elem = 3;

            double[] SH = new double[n_elem];
            double[] SF = new double[n_elem];

            string tmp_SH;      string tmp_SF;

            for(int idx = 0; idx < n_elem; idx++)
            {
                tmp_SH = this.KS.GetVar($"ZPP[{idx}].Flanke.SH");
                tmp_SF = this.KS.GetVar($"ZPP[{idx}].Fuss.SFnorm");

                SH[idx] = Convert.ToDouble(tmp_SH);
                SF[idx] = Convert.ToDouble(tmp_SF);

                Console.WriteLine(tmp_SF);
            }

            double[] ret = SH.Concat(SF).ToArray();
            return ret;
        }

        public void test_gs()
        {
            Console.WriteLine("test Gear_Set.");
        }

    }
}
