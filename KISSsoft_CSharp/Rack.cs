using System;

namespace KISSsoft_CSharp
{
    class Rack
    {
        // attributes
        protected string type { get; set; }
        protected double m_n { get; set; }
        protected double alpha_n { get; set; }

        // methods
        public Rack(string tp, double mn, double alp)
        {
            this.type = tp;
            this.m_n = mn;
            this.alpha_n = alp * Math.PI / 180.0;
        }
    }
}
