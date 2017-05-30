using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ScagnosticsSharp
{
    public class BinnedData
    {
        private Int32[] x = null;
        private Int32[] y = null;
        private Int32[] counts = null;

        private const Double RESOLUTION = 1000;

        public BinnedData(Double[] x, Double[] y, Int32[] counts)
        {
            this.x = Int32egerizeData(x);
            this.y = Int32egerizeData(y);
            this.counts = counts;
        }

        private Int32[] Int32egerizeData(Double[] x)
        {
            Int32 n = x.Length;
            Int32[] xd = new Int32[n];
            for (Int32 i = 0; i < n; i++)
            {
                xd[i] = (Int32)(RESOLUTION * x[i]);
            }
            return xd;
        }

        public Int32[] getXData()
        {
            return x;
        }

        public Int32[] getYData()
        {
            return y;
        }

        public Int32[] getCounts()
        {
            return counts;
        }
    }
}
