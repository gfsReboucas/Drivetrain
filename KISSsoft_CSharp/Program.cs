using System;

namespace KISSsoft_CSharp
{
    class Program
    {
        static void Main(string[] args)
        {
            Console.WriteLine("Hello World!");
            int[] zg = { 19 };
            double[] xg = { 0.5 };
            Gear g = new Gear(45.0, 20.0, "A", zg, 100.0, xg, 0.0);
            g.test();

            Console.ReadKey();
        }
    }
}
