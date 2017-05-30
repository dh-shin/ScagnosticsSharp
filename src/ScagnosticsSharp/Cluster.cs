using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ScagnosticsSharp
{
    public class Cluster
    {
        private Int32[] members;
        private Int32 numClusters;
        private Int32 numIterations;
        private Int32 nVar;
        private Int32 nRow;

        public Cluster(Int32 numClusters, Int32 numIterations)
        {
            this.numIterations = 3;
            this.numClusters = 0;
            if (numIterations != 0) this.numIterations = numIterations;
            if (numClusters != 0) this.numClusters = numClusters;
        }

        public Int32[] compute(Double[,] data)
        {
            nRow = data.GetLength(0);
            nVar = data.GetLength(1);
            Boolean useStoppingRule = false;
            Double[,] ssr = null;
            if (numClusters == 0)
            {
                useStoppingRule = true;
                numClusters = 25;
                ssr = new Double[numClusters, nVar];
            }

            Double[,] center = new Double[numClusters, nVar];
            Double[,] count = new Double[numClusters, nVar];
            Double[,] mean = new Double[numClusters, nVar];
            Double[,] min = new Double[numClusters, nVar];
            Double[,] max = new Double[numClusters, nVar];
            Double[,] ssq = new Double[numClusters, nVar];
            Int32[] closestPoInt32s = new Int32[numClusters];
            Double[] closestDistances = new Double[numClusters];

            members = new Int32[nRow];
            Int32[] mem = new Int32[nRow];

            for (Int32 k = 0; k < numClusters; k++)
            {

                /* find best assignments for current number of clusters */
                for (Int32 iter = 0; iter < numIterations; iter++)
                {
                    Boolean reassigned = false;
                    for (Int32 l = 0; l <= k; l++)
                    {
                        for (Int32 j = 0; j < nVar; j++)
                        {
                            if (iter == 0 || center[l, j] != mean[l, j])
                            {
                                reassign(k, data, center, count, mean, min, max, ssq,
                                        closestPoInt32s, closestDistances);
                                reassigned = true;
                                break;
                            }
                        }
                        if (reassigned)
                            break;
                    }
                    if (!reassigned || k == 0)
                        break;
                }

                /* check whether we should stop */

                if (useStoppingRule)
                {
                    Double ssq1 = 0;
                    Double ssq2 = 0;
                    for (Int32 j = 0; j < nVar; j++)
                    {
                        for (Int32 l = 0; l <= k; l++)
                        {
                            ssq1 += ssr[l, j];
                            ssq2 += ssq[l, j];
                            ssr[l, j] = ssq[l, j];
                        }
                    }
                    Double pre = (ssq1 - ssq2) / ssq1;       //proportional reduction of error
                    if (pre > 0 && pre < .1)
                    {
                        numClusters = k;
                        reassign(k, data, center, count, mean, min, max, ssq,
                                closestPoInt32s, closestDistances);
                        Array.Copy(mem, members, nRow);
                        break;
                    }
                    else
                    {
                        Array.Copy(members, mem, nRow);
                    }
                }

                /* now split a cluster to increment number of clusters */

                if (k < numClusters - 1)
                {
                    Int32 kn = k + 1;
                    Double dmax = 0;
                    Int32 jm = 0;
                    Int32 km = 0;
                    Double cutpoInt32 = 0;
                    for (Int32 l = 0; l <= k; l++)
                    {
                        for (Int32 j = 0; j < nVar; j++)
                        {
                            Double dm = max[l, j] - min[l, j];
                            if (dm > dmax)
                            {
                                cutpoInt32 = mean[l, j];
                                dmax = dm;
                                jm = j;
                                km = l;
                            }
                        }
                    }
                    for (Int32 i = 0; i < nRow; i++)
                    {
                        if (members[i] == km && data[i, jm] > cutpoInt32)
                        {
                            for (Int32 j = 0; j < nVar; j++)
                            {
                                count[km, j]--;
                                count[kn, j]++;
                                mean[km, j] -= (data[i, j] - mean[km, j]) / count[km, j];
                                mean[kn, j] += (data[i, j] - mean[kn, j]) / count[kn, j];
                            }
                        }
                    }
                }
            }

            Int32 nc = 0;
            Double cutoff = .1;
            for (Int32 k = 0; k < numClusters; k++)
            {
                if (count[k, 0] / (Double)nRow > cutoff) nc++;
            }

            Int32[] exemplars = new Int32[nc];
            nc = 0;
            for (Int32 k = 0; k < numClusters; k++)
            {
                if (count[k, 0] / (Double)nRow > cutoff)
                {
                    exemplars[nc] = closestPoInt32s[k];
                    nc++;
                }
            }
            return exemplars;
        }

        private void reassign(Int32 nCluster, Double[,] data, Double[,] center, Double[,] count,
                              Double[,] mean, Double[,] min, Double[,] max, Double[,] ssq,
                              Int32[] closestPoInt32s, Double[] closestDistances)
        {

            /* initialize cluster statistics */

            for (Int32 k = 0; k <= nCluster; k++)
            {
                closestPoInt32s[k] = -1;
                closestDistances[k] = Double.PositiveInfinity;
                for (Int32 j = 0; j < nVar; j++)
                {
                    center[k, j] = mean[k, j];
                    mean[k, j] = 0;
                    count[k, j] = 0;
                    ssq[k, j] = 0;
                    min[k, j] = Double.PositiveInfinity;
                    max[k, j] = Double.NegativeInfinity;
                }
            }

            /* assign each poInt32 to closest cluster and update statistics */


            for (Int32 i = 0; i < nRow; i++)
            {

                /* find closest cluster */
                Double dmin = Double.PositiveInfinity;
                Int32 kmin = -1;
                for (Int32 k = 0; k <= nCluster; k++)
                {
                    Double dd = distance(data, center, i, k);
                    if (dd < dmin)
                    {
                        dmin = dd;
                        kmin = k;
                        if (dmin < closestDistances[k])
                        {
                            closestDistances[k] = dmin;
                            closestPoInt32s[k] = i;
                        }
                    }
                }
                if (kmin < 0)
                {
                    members[i] = -1;
                }
                else
                {
                    members[i] = kmin;
                }

                /* update cluster statistics */

                for (Int32 j = 0; j < nVar; j++)
                {
                    if (!Double.IsNaN(data[i, j]))
                    {
                        count[kmin, j]++;
                        Double xn = count[kmin, j];
                        Double xa = data[i, j];
                        mean[kmin, j] += (xa - mean[kmin, j]) / xn;
                        if (xn > 1.0)
                            ssq[kmin, j] += xn * (xa - mean[kmin, j]) * (xa - mean[kmin, j]) / (xn - 1.0);
                        if (min[kmin, j] > xa)
                            min[kmin, j] = xa;
                        if (max[kmin, j] < xa)
                            max[kmin, j] = xa;
                    }
                }
            }
        }

        private Double distance(Double[,] a, Double[,] b, Int32 ainx, Int32 binx)
        {
            Double dist = 0;
            for (Int32 i = 0; i < a.GetLength(0); i++)
            {
                dist += (a[ainx, i] - b[binx, i]) * (a[ainx, i] - b[binx, i]);
            }
            return dist;
        }
    }
}
