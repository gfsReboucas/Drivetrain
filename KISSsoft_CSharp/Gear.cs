using System;

namespace KISSsoft_CSharp
{
    class Gear:Rack
    {
        // attributes
        protected int[] z;
        protected double b;
        protected double[] x;
        protected double beta;

        public Gear(double mn, double alp, string tp, int[] zz, double bb, double[] xx, double bet) : base(tp, mn, alp)
        {
            this.z = zz;
            this.b = bb;
            this.x = xx;
            this.beta = bet * Math.PI / 180.0;
        }

        public void test()
        {
            Console.WriteLine("test Gear");
        }
    }
}
