using System;

namespace KISSsoft_CSharp
{
    class Program
    {
        static void Main(string[] args)
        {
            int[] zg = { 19, 17, 56 };
            double[] xg = { 0.617, 0.802, -0.501 };

            double[] SHF = new double[6];
            Gear_Set S1 = new Gear_Set("planetary", 3, 863.0, 45.0, 20.0, "A", zg, 491.0, xg, 0.0);
            SHF = S1.structural_integrity(5.0e3, 12.1, 1.25, 175200.0);
            S1.test();
            S1.test_gs();

            for (int idx = 0; idx < 6; idx++)
                Console.WriteLine(SHF[idx]);

            Console.ReadKey();
        }
    }
}
