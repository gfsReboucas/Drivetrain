using System;

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

            if (this.configuration == "parallel")
            {
                this.N_p = 1;

                ks.GetModule("Z012", false);

                for (int idx = 0; idx < 2; idx++)
                {
                    ks.SetVar($"ZR[{idx}].z", this.z[idx].ToString("G5"));
                    ks.SetVar($"ZR[{idx}].x.nul", this.x[idx].ToString("G5"));
                    ks.SetVar($"ZR[{idx}].b", this.b.ToString("G5"));
                }

            }
            else if(this.configuration == "planetary")
            {
                this.N_p = np;

                ks.GetModule("Z014", false);

                for (int idx = 0; idx < 3; idx++)
                {
                    ks.SetVar($"ZR[{idx}].z", this.z[idx].ToString("G5"));
                    ks.SetVar($"ZR[{idx}].x.nul", this.x[idx].ToString("G5"));
                    ks.SetVar($"ZR[{idx}].b", this.b.ToString("G5"));
                }
            }

            ks.SetVar("ZR[0].Tool.RefProfile.name", "1.25 / 0.38 / 1.0 ISO 53:1998 Profil A");

            ks.SetVar("ZS.Geo.mn"  , this.m_n.ToString("G5"));
            ks.SetVar("ZP[0].a"    , this.a_w.ToString("G5"));
            ks.SetVar("ZS.Geo.alfn", this.alpha_n.ToString("G5"));
            ks.SetVar("ZS.Geo.beta", this.beta.ToString("G5"));

            this.KS = ks;
        }

        public double[] structural_integrity(double P_r, double n_r, double K_A)
        {
            KS.SetVar("ZS.Pnominal", P_r.ToString("G5"));
            KS.SetVar("ZS.nnominal", n_r.ToString("G5"));
            KS.SetVar("ZS.KA", K_A.ToString("G5"));

            double[] ret = { P_r, n_r };
            return ret;
        }

        public void test_gs()
        {
            Console.WriteLine("test Gear_Set.");
        }

    }
}
