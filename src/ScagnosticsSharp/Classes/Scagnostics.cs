using System;
using System.Collections.Generic;
using System.Linq;
using System.IO;
using System.Reflection;

namespace ScagnosticsSharp
{
    public class Scagnostics
    {
        private const Double FUZZ = .999;
        private const Int32 NUM_BINS = 50;      // user setting for number of bins
        private const Int32 MAX_BINS = 1000;    // user setting for maximum number of nonempty bins allowed (maxBins >= numBins*numBins)
        private const Int32 DEFAULT_RANDOM_SEED = 13579;
        private const Int32 NUM_MEASURES = 9;
        private const Int32 OUTLYING = 0, SKEWED = 1, CLUMPY = 2,
            SPARSE = 3, STRIATED = 4, CONVEX = 5,
            SKINNY = 6, STRINGY = 7, MONOTONIC = 8;

        private static readonly String[] ScagnosticsLabels = 
            {"Outlying", "Skewed", "Clumpy",
            "Sparse", "Striated", "Convex",
            "Skinny", "Stringy", "Monotonic"};

        private BinnedData BinData;
        private List<Node> Nodes;           // nodes set
        private List<Edge> Edges;           // edges set
        private List<Triangle> Triangles;   // triangles set
        private List<Edge> MstEdges;        // minimum spanning tree set
        private Edge HullStart;             // entering edge of convex hull
        private Edge ActE;
        private Int32 TotalPeeledCount;
        private Int32 TotalCount;
        private Double AlphaArea = 1;
        private Double AlphaPerimeter = 1;
        private Double HullArea = 1;
        private Double HullPerimeter = 1;
        private Double TotalOriginalMSTLengths;
        private Double TotalMSTOutlierLengths;
        private Double[] SortedOriginalMSTLengths;

        private Int32[] Px, Py, Counts;
        private Boolean[] IsOutlier;

        // For Random Number for Delaunay Triangulation
        private Int32 randomSeed;
        public Int32 RandomSeed
        {
            get { return randomSeed; }
            set { randomSeed = value; }
        }
        private static Boolean IsJavaRandReady = false;
        private static Double[] Rx;
        private static Double[] Ry;

        public Scagnostics(Double[] x, Double[] y, Int32 numBins = NUM_BINS, Int32 maxBins = MAX_BINS, Boolean doNormalize = false)
        {
            Nodes = new List<Node>();
            Edges = new List<Edge>();
            Triangles = new List<Triangle>();
            MstEdges = new List<Edge>();

            if(doNormalize == true)
            {
                NormalizeArray(x);
                NormalizeArray(y);
            }

            Binner b = new Binner(maxBins);
            BinData = b.BinHex(x, y, numBins);

            RandomSeed = DEFAULT_RANDOM_SEED;
        }

        public static Int32 GetNumScagnostics()
        {
            return NUM_MEASURES;
        }

        public static String[] GetScagnosticsLabels()
        {
            return ScagnosticsLabels;
        }

        private static void NormalizeArray(Double[] data)
        {
            Double min = Double.MaxValue;
            Double max = Double.MinValue;
            for (Int32 i = 0; i < data.Length; i++)
            {
                if (i == 0)
                    max = min = data[0];
                else if (min > data[i])
                    min = data[i];
                else if (max < data[i])
                    max = data[i];
            }
            for (Int32 i = 0; i < data.Length; i++)
                data[i] = (data[i] - min) / (max - min);
        }

        public Double[] Compute()
        {
            Px = BinData.X;
            Py = BinData.Y;
            if (Px.Length < 3)
                return null;
            Int32 xx = Px[0];
            Int32 yy = Py[0];
            Boolean isXConstant = true;
            Boolean isYConstant = true;
            for (Int32 i = 1; i < Px.Length; i++)
            {
                if (Px[i] != xx) isXConstant = false;
                if (Py[i] != yy) isYConstant = false;
            }
            if (isXConstant || isYConstant)
                return null;

            FindOutliers(BinData);

            ComputeAlphaGraph();
            ComputeTotalCount();
            ComputeAlphaArea();
            ComputeAlphaPerimeter();
            ComputeHullArea();
            ComputeHullPerimeter();

            return ComputeMeasures();
        }

        public static Boolean[] ComputeScagnosticsExemplars(Double[,] pts)
        {
            Int32 nPts = pts.GetLength(0);
            if (nPts < 2)
                return null;
            Cluster c = new Cluster(0, 0);
            Int32[] exemp = c.Compute(pts);
            Boolean[] exemplars = new Boolean[nPts];
            for (Int32 i = 0; i < exemp.Length; i++)
                exemplars[exemp[i]] = true;
            return exemplars;
        }

        public static Boolean[] ComputeScagnosticsOutliers(Double[,] pts)
        {
            // Prim's algorithm
            Int32 nPts = pts.GetLength(0);  // p*(p-1)/2 poInt32s representing pairwise scatterplots
            Int32 nVar = pts.GetLength(1);  // number of scagnostics (9)
            if (nPts < 2)
                return null;
            Int32[,] edges = new Int32[nPts - 1, 2];
            Int32[] list = new Int32[nPts];
            Int32[] degrees = new Int32[nPts];
            Double[] cost = new Double[nPts];
            Double[] lengths = new Double[nPts - 1];

            list[0] = 0;
            cost[0] = Double.PositiveInfinity;
            Int32 cheapest = 0;

            for (Int32 i = 1; i < nPts; i++)
            {
                for (Int32 j = 0; j < nVar; j++)
                {
                    Double d = pts[i, j] - pts[0, j];
                    cost[i] += d * d;
                }
                if (cost[i] < cost[cheapest])
                    cheapest = i;
            }
            for (Int32 j = 1; j < nPts; j++)
            {
                Int32 end = list[cheapest];
                Int32 jp = j - 1;
                edges[jp, 0] = cheapest;
                edges[jp, 1] = end;
                lengths[jp] = cost[cheapest];
                degrees[cheapest]++;
                degrees[end]++;
                cost[cheapest] = Double.PositiveInfinity;
                end = cheapest;

                for (Int32 i = 1; i < nPts; i++)
                {
                    if (cost[i] != Double.PositiveInfinity)
                    {
                        Double dist = 0.0;
                        for (Int32 k = 0; k < nVar; k++)
                        {
                            Double d = pts[i, k] - pts[end, k];
                            dist += d * d;
                        }
                        if (dist < cost[i])
                        {
                            list[i] = end;
                            cost[i] = dist;
                        }
                        if (cost[i] < cost[cheapest]) cheapest = i;
                    }
                }
            }
            Double cutoff = FindCutoff(lengths);
            Boolean[] outliers = new Boolean[nPts];
            for (Int32 i = 0; i < nPts; i++)
                outliers[i] = true;
            for (Int32 i = 0; i < nPts - 1; i++)
            {
                if (lengths[i] < cutoff)
                {
                    for (Int32 k = 0; k < 2; k++)
                    {
                        Int32 node = edges[i, k];
                        outliers[node] = false;
                    }
                }
            }
            return outliers;
        }

        public static void LoadJavaRandomNumber()
        {
            if (IsJavaRandReady == true)
                return;

            var assembly = Assembly.GetExecutingAssembly();
            var resourceName = "ScagnosticsSharp.Resources.java_random_seed13579_row1000.txt";
            List<List<Double>> lists;
            char[] delimChars = { ',' };
            using (Stream stream = assembly.GetManifestResourceStream(resourceName))
            using (StreamReader reader = new StreamReader(stream))
            {
                lists = new List<List<Double>>();
                while (reader.EndOfStream == false)
                {
                    String result = reader.ReadLine();
                    List<Double> values = new List<Double>();
                    String[] tokenized = result.Split(delimChars);
                    for (Int32 i = 0; i < tokenized.Length; i++)
                    {
                        Double val;
                        if (Double.TryParse(tokenized[i], out val) == true)
                        {
                            values.Add(val);
                        }
                    }
                    lists.Add(values);
                }
            }

            Int32 nVar = lists.First().Count;
            Int32 nRow = lists.Count;

            Rx = new Double[nRow];
            for (Int32 i = 0; i < nRow; i++)
            {
                Rx[i] = lists[i][0];
            }

            Ry = new Double[nRow];
            for (Int32 i = 0; i < nRow; i++)
            {
                Ry[i] = lists[i][1];
            }

            IsJavaRandReady = true;
        }

        public static void UnloadJavaRandomNumber()
        {
            IsJavaRandReady = false;
        }

        private void Clear()
        {
            Nodes.Clear();
            Edges.Clear();
            Triangles.Clear();
            MstEdges.Clear();
        }

        private void FindOutliers(BinnedData bdata)
        {
            this.Counts = bdata.Counts;
            IsOutlier = new Boolean[Px.Length];
            ComputeDT(Px, Py);
            ComputeMST();
            SortedOriginalMSTLengths = GetSortedMSTEdgeLengths();
            Double cutoff = ComputeCutoff(SortedOriginalMSTLengths);
            ComputeTotalOriginalMSTLengths();
            Boolean foundNewOutliers = ComputeMSTOutliers(cutoff);
            Double[] sortedPeeledMSTLengths;
            while (foundNewOutliers)
            {
                Clear();
                ComputeDT(Px, Py);
                ComputeMST();
                sortedPeeledMSTLengths = GetSortedMSTEdgeLengths();
                cutoff = ComputeCutoff(sortedPeeledMSTLengths);
                foundNewOutliers = ComputeMSTOutliers(cutoff);
            }
        }

        private void ComputeTotalCount()
        {
            for (Int32 i = 0; i < Counts.Length; i++)
            {
                TotalCount += Counts[i];
            }
        }

        private Double[] ComputeMeasures()
        {
            Double[] results = new Double[NUM_MEASURES];

            // Do not change order of these calls!
            results[OUTLYING] = ComputeOutlierMeasure();
            results[CLUMPY] = ComputeClusterMeasure();
            results[SKEWED] = ComputeMSTEdgeLengthSkewnessMeasure();
            results[CONVEX] = ComputeConvexityMeasure();
            results[SKINNY] = ComputeSkinnyMeasure();
            results[STRINGY] = ComputeStringyMeasure();
            results[STRIATED] = ComputeStriationMeasure();
            results[SPARSE] = ComputeSparsenessMeasure();
            results[MONOTONIC] = ComputeMonotonicityMeasure();

            return results;
        }

        private void ComputeDT(Int32[] px, Int32[] py)
        {
            TotalPeeledCount = 0;
            Random rand = new Random(RandomSeed);
            for (Int32 i = 0; i < px.Length; i++)
            {
                Double rx, ry;
                if(IsJavaRandReady == true)
                {
                    rx = Rx[i];
                    ry = Ry[i];
                }
                else
                {
                    rx = rand.NextDouble();
                    ry = rand.NextDouble();
                }

                Int32 x = px[i] + (Int32)(8 * (rx - .5)); // perturb to prevent singularities
                Int32 y = py[i] + (Int32)(8 * (ry - .5)); // perturb to prevent singularities
                Int32 count = Counts[i];
                if (!IsOutlier[i])
                {
                    Insert(x, y, count, i);
                    TotalPeeledCount += count;
                }
            }
            SetNeighbors();
            MarkHull();
        }

        private void ComputeMST()
        {
            if (Nodes.Count > 1)
            {
                List<Node> mstNodes = new List<Node>();
                Node mstNode = Nodes[0];
                UpdateMSTNodes(mstNode, mstNodes);
                Int32 count = 1;
                while (count < Nodes.Count)
                {
                    Edge addEdge = null;
                    Double wmin = Double.MaxValue;
                    Node nmin = null;

                    foreach(Node n in mstNodes)
                    {
                        mstNode = n;
                        Edge candidateEdge = mstNode.ShortestEdge(false);
                        if (candidateEdge != null)
                        {
                            Double wt = candidateEdge.Weight;
                            if (wt < wmin)
                            {
                                wmin = wt;
                                nmin = mstNode;
                                addEdge = candidateEdge;
                            }
                        }
                    }
                    if (addEdge != null)
                    {
                        Node addNode = addEdge.OtherNode(nmin);
                        UpdateMSTNodes(addNode, mstNodes);
                        UpdateMSTEdges(addEdge, MstEdges);
                    }
                    count++;
                }
            }
        }

        private static Double FindCutoff(Double[] distances)
        {
            Int32[] index = Sorts.IndexedDoubleArraySort(distances, 0, 0);
            Int32 n50 = distances.Length / 2;
            Int32 n25 = n50 / 2;
            Int32 n75 = n50 + n50 / 2;
            return distances[index[n75]] + 1.5 * (distances[index[n75]] - distances[index[n25]]);
        }

        private Boolean ComputeMSTOutliers(Double omega)
        {
            Boolean found = false;

            foreach(Node n in Nodes)
            {
                Boolean delete = true;
                foreach (Edge e in n.Neighbors)
                {
                    if (e.OnMST && e.Weight < omega)
                        delete = false;
                }
                if (delete)
                {
                    Double sumlength = 0;
                    foreach (Edge e in n.Neighbors)
                    {
                        if (e.OnMST && !e.OnOutlier)
                        {
                            sumlength += e.Weight;
                            e.OnOutlier = true;
                        }
                    }
                    TotalMSTOutlierLengths += sumlength;
                    IsOutlier[n.PointID] = true;
                    found = true;
                }
            }

            return found;
        }

        private Double ComputeCutoff(Double[] lengths)
        {
            if (lengths.Length == 0) return 0;
            Int32 n50 = lengths.Length / 2;
            Int32 n25 = n50 / 2;
            Int32 n75 = n50 + n25;
            return lengths[n75] + 1.5 * (lengths[n75] - lengths[n25]);
        }

        private Double ComputeAlphaValue()
        {
            Int32 length = SortedOriginalMSTLengths.Length;
            if (length == 0) return 100.0;
            Int32 n90 = (9 * length) / 10;
            Double alpha = SortedOriginalMSTLengths[n90];
            return Math.Min(alpha, 100.0);
        }

        private Double ComputeMSTEdgeLengthSkewnessMeasure()
        {
            if (SortedOriginalMSTLengths.Length == 0)
                return 0;
            Int32 n = SortedOriginalMSTLengths.Length;
            Int32 n50 = n / 2;
            Int32 n10 = n / 10;
            Int32 n90 = (9 * n) / 10;
            Double skewness = (SortedOriginalMSTLengths[n90] - SortedOriginalMSTLengths[n50]) /
                    (SortedOriginalMSTLengths[n90] - SortedOriginalMSTLengths[n10]);
            Double t = (Double)TotalCount / 500;
            Double correction = .7 + .3 / (1 + t * t);
            return 1 - correction * (1 - skewness);
        }

        private void UpdateMSTEdges(Edge addEdge, List<Edge> mstEdges)
        {
            mstEdges.Add(addEdge);
            addEdge.OnMST = true;
            addEdge.P1.MstDegree++;
            addEdge.P2.MstDegree++;
        }

        private void UpdateMSTNodes(Node addNode, List<Node> mstNodes)
        {
            mstNodes.Add(addNode);
            addNode.OnMST = true;
        }

        private Double[] GetSortedMSTEdgeLengths()
        {
            Double[] lengths = ComputeEdgeLengths(MstEdges, MstEdges.Count);
            Sorts.DoubleArraySort(lengths, 0, 0);
            return lengths;
        }

        private void ComputeTotalOriginalMSTLengths()
        {
            for (Int32 i = 0; i < SortedOriginalMSTLengths.Length; i++)
                TotalOriginalMSTLengths += SortedOriginalMSTLengths[i];
        }

        private Double ComputeOutlierMeasure()
        {
            return TotalMSTOutlierLengths / TotalOriginalMSTLengths;
        }

        private Double[] ComputeEdgeLengths(List<Edge> graph, Int32 n)
        {
            Double[] lengths = new Double[n];
            Int32 i = 0;
            foreach(Edge e in graph)
            {
                lengths[i] = e.Weight;
                i++;
            }
            return lengths;
        }

        private Boolean PointsInCircle(Node n, Double xc, Double yc, Double radius)
        {
            Double r = FUZZ * radius;
            foreach(Edge e in n.Neighbors)
            {
                Node no = e.OtherNode(n);
                Double dist = no.DistToNode(xc, yc);
                if (dist < r)
                    return true;
            }
            return false;
        }

        private void ComputeAlphaGraph()
        { 
            // requires initializing SEdge.onShape = false
            Boolean deleted;
            Double alpha = ComputeAlphaValue();
            do
            {
                deleted = false;
                foreach (Edge e in Edges)
                {
                    if (e.InT.OnComplex)
                    {
                        if (alpha < e.Weight / 2)
                        {
                            e.InT.OnComplex = false;
                            deleted = true;
                        }
                        else
                        {
                            if (e.InvE != null)
                                if (e.InvE.InT.OnComplex)
                                    continue;
                            if (!EdgeIsExposed(alpha, e))
                            {
                                e.InT.OnComplex = false;
                                deleted = true;
                            }
                        }
                    }
                }
            } while (deleted);
            MarkShape();
        }

        private void MarkShape()
        {
            foreach(Edge e in Edges)
            {
                e.OnShape = false;
                if (e.InT.OnComplex)
                {
                    if (e.InvE == null)
                    {
                        e.OnShape = true;
                    }
                    else if (!e.InvE.InT.OnComplex)
                        e.OnShape = true;
                }
            }
        }

        private Boolean EdgeIsExposed(Double alpha, Edge e)
        {
            Double x1 = e.P1.X;
            Double x2 = e.P2.X;
            Double y1 = e.P1.Y;
            Double y2 = e.P2.Y;
            Double xe = (x1 + x2) / 2;
            Double ye = (y1 + y2) / 2;
            Double d = Math.Sqrt(alpha * alpha - e.Weight * e.Weight / 4);
            Double xt = d * (y2 - y1) / e.Weight;
            Double yt = d * (x2 - x1) / e.Weight;
            Double xc1 = xe + xt;
            Double yc1 = ye - yt;
            Double xc2 = xe - xt;
            Double yc2 = ye + yt;
            Boolean pointsInCircle1 = PointsInCircle(e.P1, xc1, yc1, alpha) ||
                    PointsInCircle(e.P2, xc1, yc1, alpha);
            Boolean pointsInCircle2 = PointsInCircle(e.P1, xc2, yc2, alpha) ||
                    PointsInCircle(e.P2, xc2, yc2, alpha);
            return !(pointsInCircle1 && pointsInCircle2);
        }

        private Double ComputeStringyMeasure()
        {
            Int32 count1 = 0;
            Int32 count2 = 0;

            foreach(Node n in Nodes)
            {
                if (n.MstDegree == 1)
                    count1++;
                if (n.MstDegree == 2)
                    count2++;
            }
            Double result = (Double)count2 / (Double)(Nodes.Count - count1);
            return result * result * result;
        }

        private Double ComputeClusterMeasure()
        {
            Double[] maxLength = new Double[1];
            Double maxValue = 0;
            foreach(Edge e in MstEdges)
            {
                ClearVisits();
                e.OnMST = false;  // break MST at this edge
                Int32 runts = e.GetRunts(maxLength);
                e.OnMST = true;   // restore this edge to MST
                if (maxLength[0] > 0)
                {
                    Double value = runts * (1 - maxLength[0] / e.Weight);
                    if (value > maxValue)
                        maxValue = value;
                }
            }
            return 2 * maxValue / TotalPeeledCount;
        }

        private void ClearVisits()
        {
            foreach(Node n in Nodes)
            {
                n.IsVisited = false;
            }
        }

        private Double ComputeMonotonicityMeasure()
        {
            Int32 n = Counts.Length;
            Double[] ax = new Double[n];
            Double[] ay = new Double[n];
            Double[] weights = new Double[n];
            for (Int32 i = 0; i < n; i++)
            {
                ax[i] = Px[i];
                ay[i] = Py[i];
                weights[i] = Counts[i];
            }
            Double[] rx = Sorts.Rank(ax);
            Double[] ry = Sorts.Rank(ay);
            Double s = ComputePearson(rx, ry, weights);
            return s * s;
        }

        private Double ComputePearson(Double[] x, Double[] y, Double[] weights)
        {
            Int32 n = x.Length;
            Double xmean = 0;
            Double ymean = 0;
            Double xx = 0;
            Double yy = 0;
            Double xy = 0;
            Double sumwt = 0;
            for (Int32 i = 0; i < n; i++)
            {
                Double wt = weights[i];
                if (wt > 0 && !IsOutlier[i])
                {
                    sumwt += wt;
                    xx += (x[i] - xmean) * wt * (x[i] - xmean);
                    yy += (y[i] - ymean) * wt * (y[i] - ymean);
                    xy += (x[i] - xmean) * wt * (y[i] - ymean);
                    xmean += (x[i] - xmean) * wt / sumwt;
                    ymean += (y[i] - ymean) * wt / sumwt;
                }
            }
            xy = xy / Math.Sqrt(xx * yy);
            return xy;
        }

        private Double ComputeSparsenessMeasure()
        {
            Int32 n = SortedOriginalMSTLengths.Length;
            Int32 n90 = (9 * n) / 10;
            Double sparse = Math.Min(SortedOriginalMSTLengths[n90] / 1000, 1);
            Double t = (Double)TotalCount / 500;
            Double correction = .7 + .3 / (1 + t * t);
            return correction * sparse;
        }

        private Double ComputeStriationMeasure()
        {
            Double numEdges = 0;

            foreach(Edge e in MstEdges)
            {
                Node n1 = e.P1;
                Node n2 = e.P2;
                if (n1.MstDegree == 2 && n2.MstDegree == 2)
                {
                    Edge e1 = GetAdjacentMSTEdge(n1, e);
                    Edge e2 = GetAdjacentMSTEdge(n2, e);
                    if (CosineOfAdjacentEdges(e, e1, n1) < -.7 && CosineOfAdjacentEdges(e, e2, n2) < -.7)
                        numEdges++;
                }
            }
            return numEdges / (Double)MstEdges.Count;
        }

        private Edge GetAdjacentMSTEdge(Node n, Edge e)
        {
            foreach(Edge et in n.Neighbors)
            {
                if (et.OnMST && e != et)
                {
                    return et;
                }
            }
            return null;
        }

        private Double CosineOfAdjacentEdges(Edge e1, Edge e2, Node n)
        {
            Double v1x = e1.OtherNode(n).X - n.X;
            Double v1y = e1.OtherNode(n).Y - n.Y;
            Double v2x = e2.OtherNode(n).X - n.X;
            Double v2y = e2.OtherNode(n).Y - n.Y;
            Double v1 = Math.Sqrt(v1x * v1x + v1y * v1y);
            Double v2 = Math.Sqrt(v2x * v2x + v2y * v2y);
            v1x = v1x / v1;
            v1y = v1y / v1;
            v2x = v2x / v2;
            v2y = v2y / v2;
            return v1x * v2x + v1y * v2y;
        }

        private Double ComputeConvexityMeasure()
        {
            if (HullArea == 0) // poInt32s in general position
                return 1;
            else
            {
                Double t = (Double)TotalCount / 500;
                Double correction = .7 + .3 / (1 + t * t);
                Double convexity = AlphaArea / HullArea;
                return correction * convexity;
            }
        }

        private Double ComputeSkinnyMeasure()
        {
            if (AlphaPerimeter > 0)
                return 1 - Math.Sqrt(4 * Math.PI * AlphaArea) / AlphaPerimeter;
            else
                return 1;
        }

        private void ComputeAlphaArea()
        {
            Double area = 0;       
            foreach(Triangle t in Triangles)
            {
                if (t.OnComplex)
                {
                    Node p1 = t.AnEdge.P1;
                    Node p2 = t.AnEdge.P2;
                    Node p3 = t.AnEdge.NextE.P2;
                    area += Math.Abs(p1.X * p2.Y + p1.Y * p3.X + p2.X * p3.Y
                            - p3.X * p2.Y - p3.Y * p1.X - p1.Y * p2.X);
                }
            }
            AlphaArea = area / 2;
        }

        private void ComputeHullArea()
        {
            Double area = 0.0;
            foreach(Triangle t in Triangles)
            {
                Node p1 = t.AnEdge.P1;
                Node p2 = t.AnEdge.P2;
                Node p3 = t.AnEdge.NextE.P2;
                area += Math.Abs(p1.X * p2.Y + p1.Y * p3.X + p2.X * p3.Y
                        - p3.X * p2.Y - p3.Y * p1.X - p1.Y * p2.X);
            }
            HullArea = area / 2.0;
        }

        private void ComputeAlphaPerimeter()
        {
            Double sum = 0;
            foreach(Edge e in Edges)
            {
                if (e.OnShape)
                {
                    sum += e.Weight;
                }
            }
            AlphaPerimeter = sum;
        }

        private void ComputeHullPerimeter()
        {
            Double sum = 0;
            Edge e = HullStart;
            do
            {
                sum += e.P1.DistToNode(e.P2.X, e.P2.Y);
                e = e.NextH;
            } while (!e.IsEqual(HullStart));
            HullPerimeter = sum;
        }

        private void SetNeighbors()
        {
            foreach(Edge e in Edges)
            {
                if (e.IsNewEdge(e.P1))
                    e.P1.SetNeighbor(e);
                if (e.IsNewEdge(e.P2))
                    e.P2.SetNeighbor(e);
            }
        }

        private void Insert(Int32 px, Int32 py, Int32 count, Int32 id)
        {
            Int32 eid;
            Node nd = new Node(px, py, count, id);
            Nodes.Add(nd);
            if (Nodes.Count < 3) return;
            if (Nodes.Count == 3)    // create the first triangle
            {
                Node p1 = Nodes[0];
                Node p2 = Nodes[1];
                Node p3 = Nodes[2];
                Edge e1 = new Edge(p1, p2);
                if (e1.OnSide(p3) == 0)
                {
                    Nodes.Remove(nd);
                    return;
                }
                if (e1.OnSide(p3) == -1)  // right side
                {
                    p1 = Nodes[1];
                    p2 = Nodes[0];
                    e1.Update(p1, p2);
                }
                Edge e2 = new Edge(p2, p3);
                Edge e3 = new Edge(p3, p1);
                e1.NextH = e2;
                e2.NextH = e3;
                e3.NextH = e1;
                HullStart = e1;
                Triangles.Add(new Triangle(Edges, e1, e2, e3));
                return;
            }
            ActE = Edges[0];
            if (ActE.OnSide(nd) == -1)
            {
                if (ActE.InvE == null)
                    eid = -1;
                else
                    eid = SearchEdge(ActE.InvE, nd);
            }
            else
                eid = SearchEdge(ActE, nd);
            if (eid == 0)
            {
                Nodes.Remove(nd);
                return;
            }
            if (eid > 0)
                ExpandTri(ActE, nd, eid);   // nd is inside or on a triangle
            else
                ExpandHull(nd);                // nd is outside convex hull
        }

        private void ExpandTri(Edge e, Node nd, Int32 type)
        {
            Edge e1 = e;
            Edge e2 = e1.NextE;
            Edge e3 = e2.NextE;
            Node p1 = e1.P1;
            Node p2 = e2.P1;
            Node p3 = e3.P1;
            if (type == 2)
            {   // nd is inside of the triangle
                Edge e10 = new Edge(p1, nd);
                Edge e20 = new Edge(p2, nd);
                Edge e30 = new Edge(p3, nd);
                e.InT.RemoveEdges(Edges);
                Triangles.Remove(e.InT);     // remove old triangle
                Edge e100 = e10.MakeSymmEdge();
                Edge e200 = e20.MakeSymmEdge();
                Edge e300 = e30.MakeSymmEdge();
                Triangles.Add(new Triangle(Edges, e1, e20, e100));
                Triangles.Add(new Triangle(Edges, e2, e30, e200));
                Triangles.Add(new Triangle(Edges, e3, e10, e300));
                SwapTest(e1);   // swap test for the three new triangles
                SwapTest(e2);
                SwapTest(e3);
            }
            else
            {   
                // nd is on the edge e
                Edge e4 = e1.InvE;
                if (e4 == null || e4.InT == null)
                {          
                    // one triangle involved
                    Edge e30 = new Edge(p3, nd);
                    Edge e02 = new Edge(nd, p2);
                    Edge e10 = new Edge(p1, nd);
                    Edge e03 = e30.MakeSymmEdge();
                    //shareEdges(e03,e30);
                    e10.AsIndex();
                    e1.MostLeft().NextH = e10;
                    e10.NextH = e02;
                    e02.NextH = e1.NextH;
                    HullStart = e02;
                    Triangles.Remove(e1.InT);  // remove oldtriangle and add two new triangles
                    Edges.Remove(e1);
                    Edges.Add(e10);
                    Edges.Add(e02);
                    Edges.Add(e30);
                    Edges.Add(e03);
                    Triangles.Add(new Triangle(e2, e30, e02));
                    Triangles.Add(new Triangle(e3, e10, e03));
                    SwapTest(e2);   // swap test for the two new triangles
                    SwapTest(e3);
                    SwapTest(e30);
                }
                else
                {        
                    // two triangle involved
                    Edge e5 = e4.NextE;
                    Edge e6 = e5.NextE;
                    Node p4 = e6.P1;
                    Edge e10 = new Edge(p1, nd);
                    Edge e20 = new Edge(p2, nd);
                    Edge e30 = new Edge(p3, nd);
                    Edge e40 = new Edge(p4, nd);
                    Triangles.Remove(e.InT);                   // remove oldtriangle
                    e.InT.RemoveEdges(Edges);
                    Triangles.Remove(e4.InT);               // remove old triangle
                    e4.InT.RemoveEdges(Edges);
                    e5.AsIndex();   // because e, e4 removed, reset edge sortOrder of node p1 and p2
                    e2.AsIndex();
                    Triangles.Add(new Triangle(Edges, e2, e30, e20.MakeSymmEdge()));
                    Triangles.Add(new Triangle(Edges, e3, e10, e30.MakeSymmEdge()));
                    Triangles.Add(new Triangle(Edges, e5, e40, e10.MakeSymmEdge()));
                    Triangles.Add(new Triangle(Edges, e6, e20, e40.MakeSymmEdge()));
                    SwapTest(e2);   // swap test for the three new triangles
                    SwapTest(e3);
                    SwapTest(e5);
                    SwapTest(e6);
                    SwapTest(e10);
                    SwapTest(e20);
                    SwapTest(e30);
                    SwapTest(e40);
                }
            }
        }

        private void ExpandHull(Node nd)
        {
            Edge e1, e2, e3 = null, enext;
            Edge e = HullStart;
            Edge comedge = null, lastbe = null;
            while (true)
            {
                enext = e.NextH;
                if (e.OnSide(nd) == -1)
                {  
                    // right side
                    if (lastbe != null)
                    {
                        e1 = e.MakeSymmEdge();
                        e2 = new Edge(e.P1, nd);
                        e3 = new Edge(nd, e.P2);
                        if (comedge == null)
                        {
                            HullStart = lastbe;
                            lastbe.NextH = e2;
                            lastbe = e2;
                        }
                        else
                            comedge.LinkSymmEdge(e2);

                        comedge = e3;
                        Triangles.Add(new Triangle(Edges, e1, e2, e3));
                        SwapTest(e);
                    }
                }
                else
                {
                    if (comedge != null) break;
                    lastbe = e;
                }
                e = enext;
            }

            lastbe.NextH = e3;
            e3.NextH = e;
        }

        private Int32 SearchEdge(Edge e, Node nd)
        {
            Int32 f2, f3;
            Edge e0 = null;
            if ((f2 = e.NextE.OnSide(nd)) == -1)
            {
                if (e.NextE.InvE != null)
                    return SearchEdge(e.NextE.InvE, nd);
                else
                {
                    ActE = e;
                    return -1;
                }
            }
            if (f2 == 0) e0 = e.NextE;
            Edge ee = e.NextE;
            if ((f3 = ee.NextE.OnSide(nd)) == -1)
            {
                if (ee.NextE.InvE != null)
                    return SearchEdge(ee.NextE.InvE, nd);
                else
                {
                    ActE = ee.NextE;
                    return -1;
                }
            }
            if (f3 == 0) e0 = ee.NextE;
            if (e.OnSide(nd) == 0) e0 = e;
            if (e0 != null)
            {
                ActE = e0;
                if (e0.NextE.OnSide(nd) == 0)
                {
                    ActE = e0.NextE;
                    return 0;
                }
                if (e0.NextE.NextE.OnSide(nd) == 0) return 0;
                return 1;
            }
            ActE = ee;
            return 2;
        }

        private void SwapTest(Edge e11)
        {
            Edge e21 = e11.InvE;
            if (e21 == null || e21.InT == null) return;
            Edge e12 = e11.NextE;
            Edge e13 = e12.NextE;
            Edge e22 = e21.NextE;
            Edge e23 = e22.NextE;
            if (e11.InT.InCircle(e22.P2) || e21.InT.InCircle(e12.P2))
            {
                e11.Update(e22.P2, e12.P2);
                e21.Update(e12.P2, e22.P2);
                e11.LinkSymmEdge(e21);
                e13.InT.Update(e13, e22, e11);
                e23.InT.Update(e23, e12, e21);
                e12.AsIndex();
                e22.AsIndex();
                SwapTest(e12);
                SwapTest(e22);
                SwapTest(e13);
                SwapTest(e23);
            }
        }

        private void MarkHull()
        {
            Edge e = HullStart;
            if (e != null)
                do
                {
                    e.OnHull = true;
                    e.P1.OnHull = true;
                    e.P2.OnHull = true;
                    e = e.NextH;
                } while (!e.IsEqual(HullStart));
        }
    }
}
