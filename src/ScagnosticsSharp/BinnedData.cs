using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ScagnosticsSharp
{
    internal class BinnedData
    {
        private const Double RESOLUTION = 1000;

        public Int32[] X = null;
        public Int32[] Y = null;
        public Int32[] Counts = null;

        public BinnedData(Double[] x, Double[] y, Int32[] counts)
        {
            X = IntegerizeData(x);
            Y = IntegerizeData(y);
            Counts = counts;
        }

        private Int32[] IntegerizeData(Double[] x)
        {
            Int32 n = x.Length;
            Int32[] xd = new Int32[n];
            for (Int32 i = 0; i < n; i++)
            {
                xd[i] = (Int32)(RESOLUTION * x[i]);
            }
            return xd;
        }
    }
}
