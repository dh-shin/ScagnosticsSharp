using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ScagnosticsSharp
{
    public class Binner
    {
        private Int32 maxBins;

        public Binner(Int32 maxBins)
        {
            this.maxBins = maxBins;
        }

        public BinnedData binHex(Double[] x, Double[] y, Int32 nBins)
        {
            Int32 n = x.Length;

            // scaling constants

            Double con1 = .25;
            Double con2 = 1.0 / 3.0;
            Double c1 = (Double)(nBins - 1);
            Double c2 = c1 / Math.Sqrt(3.0);
            Int32 jinc = nBins;
            Int32 iinc = 2 * nBins;
            Int32 nBin = (nBins + 20) * (nBins + 20);

            Int32[] count = new Int32[nBin];
            Double[] xbin = new Double[nBin];
            Double[] ybin = new Double[nBin];

            // fill bins

            for (Int32 i = 0; i < n; i++)
            {
                if (Double.IsNaN(x[i])) continue;
                if (Double.IsNaN(y[i])) continue;
                Double sx = c1 * x[i];
                Double sy = c2 * y[i];
                Int32 i1 = (Int32)(sy + .5);
                Int32 j1 = (Int32)(sx + .5);
                Double dy = sy - (Double)i1;
                Double dx = sx - (Double)j1;
                Double dist1 = dx * dx + 3.0 * dy * dy;
                Int32 m;
                if (dist1 < con1)
                {
                    m = i1 * iinc + j1;
                }
                else if (dist1 > con2)
                {
                    m = (Int32)sy * iinc + (Int32)sx + jinc;
                }
                else
                {
                    Int32 i2 = (Int32)sy;
                    Int32 j2 = (Int32)sx;
                    dy = sy - (Double)i2 - .5;
                    dx = sx - (Double)j2 - .5;
                    Double dist2 = dx * dx + 3.0 * dy * dy;
                    if (dist1 <= dist2)
                    {
                        m = i1 * iinc + j1;
                    }
                    else
                    {
                        m = i2 * iinc + j2 + jinc;
                    }
                }
                count[m]++;
                xbin[m] += (x[i] - xbin[m]) / count[m];
                ybin[m] += (y[i] - ybin[m]) / count[m];
            }

            nBin = deleteEmptyBins(count, xbin, ybin);
            if (nBin > maxBins)
            {
                nBins = 2 * nBins / 3;
                return binHex(x, y, nBins);
            }

            Int32[] tcount = new Int32[nBin];
            Double[] xtbin = new Double[nBin];
            Double[] ytbin = new Double[nBin];

            Array.Copy(count, tcount, nBin);
            Array.Copy(xbin, xtbin, nBin);
            Array.Copy(ybin, ytbin, nBin);

            return new BinnedData(xtbin, ytbin, tcount);
        }

        private Int32 deleteEmptyBins(Int32[] count, Double[] xbin, Double[] ybin)
        {

            Int32 k = 0;
            for (Int32 i = 0; i < count.Length; i++)
            {
                if (count[i] > 0)
                {
                    count[k] = count[i];
                    xbin[k] = xbin[i];
                    ybin[k] = ybin[i];
                    k++;
                }
            }
            return k;
        }
    }
}
