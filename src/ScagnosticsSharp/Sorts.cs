using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ScagnosticsSharp
{
    public class Sorts
    {
        private Sorts()
        {
        }

        public static void DoubleArraySort(Double[] x, Int32 fromIndex, Int32 toIndex)
        {
            if (fromIndex == toIndex)
            {
                fromIndex = 0;
                toIndex = x.Length;
            }
            Array.Sort(x, fromIndex, toIndex);
        }

        public static Int32[] indexedDoubleArraySort(Double[] x, Int32 fromIndex, Int32 toIndex)
        {
            if (fromIndex == toIndex)
            {
                fromIndex = 0;
                toIndex = x.Length;
            }

            Double[] x_copy = new Double[x.Length];
            Array.Copy(x, x_copy, x.Length);

            Int32[] sortOrder = new Int32[toIndex - fromIndex];
            for(Int32 i = 0; i < sortOrder.Length; i++)
                sortOrder[i] = i;

            Array.Sort(x_copy, sortOrder, fromIndex, toIndex);

            return sortOrder;
        }

        public static Double[] rank(Double[] a)
        {
            Int32 k, k1, k2, kind, kms, l, lind, n;
            Double ak, am, freq;
            Boolean insert;

            n = a.Length;
            Double[] ranks = new Double[n];
            Int32[] index = indexedDoubleArraySort(a, 0, n);

            lind = index[0];
            am = a[lind];
            k1 = 0;
            k2 = 0;
            ak = 1.0;
            /* kms allows for missing data */
            kms = 1;
            for (k = 1; k < n; k++)
            {
                kind = index[k];
                insert = true;
                if (!Double.IsNaN(am))
                {
                    freq = 1.0;
                    /*
                                    if (wt != null)
                                        freq = Math.floor(wt[kind]);
                    */
                    kms += (Int32)freq;
                    if (freq > 1.0)
                    {
                        ak += 0.5 * (freq - 1.0);
                        k1 = k;
                        k2 = k;
                    }
                    else if (a[kind] == am)
                    {
                        k2 = k;
                        ak += 0.5;
                        if (k < n - 1)
                            insert = false;
                    }
                    if (insert)
                    {
                        for (l = k1; l <= k2; l++)
                        {
                            lind = index[l];
                            ranks[lind] = ak;
                        }
                        if (k2 != n - 1 && k == n - 1)
                            ranks[kind] = kms;
                    }
                }
                if (insert)
                {
                    k1 = k;
                    k2 = k;
                    ak = kms;
                    am = a[kind];
                }
            }
            return ranks;
        }
    }
}
