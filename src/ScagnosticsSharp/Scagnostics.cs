using System;
using System.Collections.Generic;
using System.Linq;
using System.IO;
using System.Reflection;

namespace ScagnosticsSharp
{
    public class Scagnostics
    {
        private BinnedData bdata;
        private List<Node> nodes;        // nodes set
        private List<Edge> edges;        // edges set
        private List<Triangle> triangles;    // triangles set
        private List<Edge> mstEdges;     // minimum spanning tree set
        private Edge hullStart;   // entering edge of convex hull
        private Edge actE;
        private Int32 totalPeeledCount;
        private Int32 totalCount;
        private Double alphaArea = 1, alphaPerimeter = 1, hullArea = 1, hullPerimeter = 1;
        private Double totalOriginalMSTLengths;
        private Double totalMSTOutlierLengths;
        private Double[] sortedOriginalMSTLengths;
        private static Int32 numScagnostics = 9;
        private static readonly Int32 OUTLYING = 0, SKEWED = 1, CLUMPY = 2, SPARSE = 3,
                STRIATED = 4, CONVEX = 5, SKINNY = 6, STRINGY = 7, MONOTONIC = 8;
        private static readonly String[] scagnosticsLabels = {"Outlying", "Skewed", "Clumpy", "Sparse",
            "Striated", "Convex", "Skinny", "Stringy", "Monotonic"};
        private Int32[] px, py, counts;
        private Boolean[] isOutlier;
        private Double FUZZ = .999;

        // For Java Random Values
        private static Boolean IsJavaRandReady = false;
        private static Double[] Rx;
        private static Double[] Ry;

        public Scagnostics(Double[] x, Double[] y, Int32 numBins, Int32 maxBins)
        {
            nodes = new List<Node>();
            edges = new List<Edge>();
            triangles = new List<Triangle>();
            mstEdges = new List<Edge>();
            Binner b = new Binner(maxBins);
            bdata = b.binHex(x, y, numBins);
        }

        public Double[] compute()
        {
            px = bdata.getXData();
            py = bdata.getYData();
            if (px.Length < 3)
                return null;
            Int32 xx = px[0];
            Int32 yy = py[0];
            Boolean isXConstant = true;
            Boolean isYConstant = true;
            for (Int32 i = 1; i < px.Length; i++)
            {
                if (px[i] != xx) isXConstant = false;
                if (py[i] != yy) isYConstant = false;
            }
            if (isXConstant || isYConstant)
                return null;

            findOutliers(bdata);

            computeAlphaGraph();
            computeTotalCount();
            computeAlphaArea();
            computeAlphaPerimeter();
            computeHullArea();
            computeHullPerimeter();
            return computeMeasures();
        }

        public static Int32 getNumScagnostics()
        {
            return scagnosticsLabels.Length;
        }

        public static String[] getScagnosticsLabels()
        {
            return scagnosticsLabels;
        }

        public static Boolean[] computeScagnosticsExemplars(Double[,] pts)
        {
            Int32 nPts = pts.GetLength(0);
            if (nPts < 2)
                return null;
            Cluster c = new Cluster(0, 0);
            Int32[] exemp = c.compute(pts);
            Boolean[] exemplars = new Boolean[nPts];
            for (Int32 i = 0; i < exemp.Length; i++)
                exemplars[exemp[i]] = true;
            return exemplars;
        }

        public static Boolean[] computeScagnosticsOutliers(Double[,] pts)
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
            Double cutoff = findCutoff(lengths);
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

        public static void LoadJavaRand()
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

        private void clear()
        {
            nodes.Clear();
            edges.Clear();
            triangles.Clear();
            mstEdges.Clear();
        }

        private void findOutliers(BinnedData bdata)
        {
            this.counts = bdata.getCounts();
            isOutlier = new Boolean[px.Length];
            computeDT(px, py);
            computeMST();
            sortedOriginalMSTLengths = getSortedMSTEdgeLengths();
            Double cutoff = computeCutoff(sortedOriginalMSTLengths);
            computeTotalOriginalMSTLengths();
            Boolean foundNewOutliers = computeMSTOutliers(cutoff);
            Double[] sortedPeeledMSTLengths;
            while (foundNewOutliers)
            {
                clear();
                computeDT(px, py);
                computeMST();
                sortedPeeledMSTLengths = getSortedMSTEdgeLengths();
                cutoff = computeCutoff(sortedPeeledMSTLengths);
                foundNewOutliers = computeMSTOutliers(cutoff);
            }
        }

        private void computeTotalCount()
        {
            for (Int32 i = 0; i < counts.Length; i++)
            {
                totalCount += counts[i];
            }
        }

        private Double[] computeMeasures()
        {
            Double[] results = new Double[numScagnostics];

            // Do not change order of these calls!
            results[OUTLYING] = computeOutlierMeasure();
            results[CLUMPY] = computeClusterMeasure();
            results[SKEWED] = computeMSTEdgeLengthSkewnessMeasure();
            results[CONVEX] = computeConvexityMeasure();
            results[SKINNY] = computeSkinnyMeasure();
            results[STRINGY] = computeStringyMeasure();
            results[STRIATED] = computeStriationMeasure();
            results[SPARSE] = computeSparsenessMeasure();
            results[MONOTONIC] = computeMonotonicityMeasure();
            return results;
        }

        private void computeDT(Int32[] px, Int32[] py)
        {
            totalPeeledCount = 0;
            Random r = new Random(13579);
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
                    rx = r.NextDouble();
                    ry = r.NextDouble();
                }

                Int32 x = px[i] + (Int32)(8 * (rx - .5)); // perturb to prevent singularities
                Int32 y = py[i] + (Int32)(8 * (ry - .5));
                Int32 count = counts[i];
                if (!isOutlier[i])
                {
                    insert(x, y, count, i);
                    totalPeeledCount += count;
                }
            }
            setNeighbors();
            markHull();
        }

        private void computeMST()
        {
            if (nodes.Count > 1)
            {
                List<Node> mstNodes = new List<Node>();
                Node mstNode = nodes[0];
                updateMSTNodes(mstNode, mstNodes);
                Int32 count = 1;
                while (count < nodes.Count)
                {
                    Edge addEdge = null;
                    Double wmin = Double.MaxValue;
                    Node nmin = null;

                    foreach(Node n in mstNodes)
                    {
                        mstNode = n;
                        Edge candidateEdge = mstNode.shortestEdge(false);
                        if (candidateEdge != null)
                        {
                            Double wt = candidateEdge.weight;
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
                        Node addNode = addEdge.otherNode(nmin);
                        updateMSTNodes(addNode, mstNodes);
                        updateMSTEdges(addEdge, mstEdges);
                    }
                    count++;
                }
            }
        }

        private static Double findCutoff(Double[] distances)
        {
            Int32[] index = Sorts.indexedDoubleArraySort(distances, 0, 0);
            Int32 n50 = distances.Length / 2;
            Int32 n25 = n50 / 2;
            Int32 n75 = n50 + n50 / 2;
            return distances[index[n75]] + 1.5 * (distances[index[n75]] - distances[index[n25]]);
        }

        private Boolean computeMSTOutliers(Double omega)
        {
            Boolean found = false;

            foreach(Node n in nodes)
            {
                Boolean delete = true;
                foreach (Edge e in n.neighbors)
                {
                    if (e.onMST && e.weight < omega)
                        delete = false;
                }
                if (delete)
                {
                    Double sumlength = 0;
                    foreach (Edge e in n.neighbors)
                    {
                        if (e.onMST && !e.onOutlier)
                        {
                            sumlength += e.weight;
                            e.onOutlier = true;
                        }
                    }
                    totalMSTOutlierLengths += sumlength;
                    isOutlier[n.poInt32ID] = true;
                    found = true;
                }
            }

            return found;
        }

        private Double computeCutoff(Double[] lengths)
        {
            if (lengths.Length == 0) return 0;
            Int32 n50 = lengths.Length / 2;
            Int32 n25 = n50 / 2;
            Int32 n75 = n50 + n25;
            return lengths[n75] + 1.5 * (lengths[n75] - lengths[n25]);
        }

        private Double computeAlphaValue()
        {
            Int32 length = sortedOriginalMSTLengths.Length;
            if (length == 0) return 100.0;
            Int32 n90 = (9 * length) / 10;
            Double alpha = sortedOriginalMSTLengths[n90];
            return Math.Min(alpha, 100.0);
        }

        private Double computeMSTEdgeLengthSkewnessMeasure()
        {
            if (sortedOriginalMSTLengths.Length == 0)
                return 0;
            Int32 n = sortedOriginalMSTLengths.Length;
            Int32 n50 = n / 2;
            Int32 n10 = n / 10;
            Int32 n90 = (9 * n) / 10;
            Double skewness = (sortedOriginalMSTLengths[n90] - sortedOriginalMSTLengths[n50]) /
                    (sortedOriginalMSTLengths[n90] - sortedOriginalMSTLengths[n10]);
            Double t = (Double)totalCount / 500;
            Double correction = .7 + .3 / (1 + t * t);
            return 1 - correction * (1 - skewness);
        }

        private void updateMSTEdges(Edge addEdge, List<Edge> mstEdges)
        {
            mstEdges.Add(addEdge);
            addEdge.onMST = true;
            addEdge.p1.mstDegree++;
            addEdge.p2.mstDegree++;
        }

        private void updateMSTNodes(Node addNode, List<Node> mstNodes)
        {
            mstNodes.Add(addNode);
            addNode.onMST = true;
        }

        private Double[] getSortedMSTEdgeLengths()
        {
            Double[] lengths = computeEdgeLengths(mstEdges, mstEdges.Count);
            Sorts.DoubleArraySort(lengths, 0, 0);
            return lengths;
        }

        private void computeTotalOriginalMSTLengths()
        {
            for (Int32 i = 0; i < sortedOriginalMSTLengths.Length; i++)
                totalOriginalMSTLengths += sortedOriginalMSTLengths[i];
        }

        private Double computeOutlierMeasure()
        {
            return totalMSTOutlierLengths / totalOriginalMSTLengths;
        }

        private Double[] computeEdgeLengths(List<Edge> graph, Int32 n)
        {
            Double[] lengths = new Double[n];
            Int32 i = 0;
            foreach(Edge e in graph)
            {
                lengths[i] = e.weight;
                i++;
            }
            return lengths;
        }

        private Boolean poInt32sInCircle(Node n, Double xc, Double yc, Double radius)
        {
            Double r = FUZZ * radius;
            foreach(Edge e in n.neighbors)
            {
                Node no = e.otherNode(n);
                Double dist = no.distToNode(xc, yc);
                if (dist < r)
                    return true;
            }
            return false;
        }

        private void computeAlphaGraph()
        { // requires initializing SEdge.onShape = false
            Boolean deleted;
            Double alpha = computeAlphaValue();
            do
            {
                deleted = false;
                foreach (Edge e in edges)
                {
                    if (e.inT.onComplex)
                    {
                        if (alpha < e.weight / 2)
                        {
                            e.inT.onComplex = false;
                            deleted = true;
                        }
                        else
                        {
                            if (e.invE != null)
                                if (e.invE.inT.onComplex)
                                    continue;
                            if (!edgeIsExposed(alpha, e))
                            {
                                e.inT.onComplex = false;
                                deleted = true;
                            }
                        }
                    }
                }
            } while (deleted);
            markShape();
        }

        private void markShape()
        {
            foreach(Edge e in edges)
            {
                e.onShape = false;
                if (e.inT.onComplex)
                {
                    if (e.invE == null)
                    {
                        e.onShape = true;
                    }
                    else if (!e.invE.inT.onComplex)
                        e.onShape = true;
                }
            }
        }

        private Boolean edgeIsExposed(Double alpha, Edge e)
        {
            Double x1 = e.p1.x;
            Double x2 = e.p2.x;
            Double y1 = e.p1.y;
            Double y2 = e.p2.y;
            Double xe = (x1 + x2) / 2;
            Double ye = (y1 + y2) / 2;
            Double d = Math.Sqrt(alpha * alpha - e.weight * e.weight / 4);
            Double xt = d * (y2 - y1) / e.weight;
            Double yt = d * (x2 - x1) / e.weight;
            Double xc1 = xe + xt;
            Double yc1 = ye - yt;
            Double xc2 = xe - xt;
            Double yc2 = ye + yt;
            Boolean poInt32sInCircle1 = poInt32sInCircle(e.p1, xc1, yc1, alpha) ||
                    poInt32sInCircle(e.p2, xc1, yc1, alpha);
            Boolean poInt32sInCircle2 = poInt32sInCircle(e.p1, xc2, yc2, alpha) ||
                    poInt32sInCircle(e.p2, xc2, yc2, alpha);
            return !(poInt32sInCircle1 && poInt32sInCircle2);
        }

        private Double computeStringyMeasure()
        {
            Int32 count1 = 0;
            Int32 count2 = 0;

            foreach(Node n in nodes)
            {
                if (n.mstDegree == 1)
                    count1++;
                if (n.mstDegree == 2)
                    count2++;
            }
            Double result = (Double)count2 / (Double)(nodes.Count - count1);
            return result * result * result;
        }

        private Double computeClusterMeasure()
        {
            Double[] maxLength = new Double[1];
            Double maxValue = 0;
            foreach(Edge e in mstEdges)
            {
                clearVisits();
                e.onMST = false;  // break MST at this edge
                Int32 runts = e.getRunts(maxLength);
                e.onMST = true;   // restore this edge to MST
                if (maxLength[0] > 0)
                {
                    Double value = runts * (1 - maxLength[0] / e.weight);
                    if (value > maxValue)
                        maxValue = value;
                }
            }
            return 2 * maxValue / totalPeeledCount;
        }

        private void clearVisits()
        {
            foreach(Node n in nodes)
            {
                n.isVisited = false;
            }
        }

        private Double computeMonotonicityMeasure()
        {
            Int32 n = counts.Length;
            Double[] ax = new Double[n];
            Double[] ay = new Double[n];
            Double[] weights = new Double[n];
            for (Int32 i = 0; i < n; i++)
            {
                ax[i] = px[i];
                ay[i] = py[i];
                weights[i] = counts[i];
            }
            Double[] rx = Sorts.rank(ax);
            Double[] ry = Sorts.rank(ay);
            Double s = computePearson(rx, ry, weights);
            return s * s;
        }

        private Double computePearson(Double[] x, Double[] y, Double[] weights)
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
                if (wt > 0 && !isOutlier[i])
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

        private Double computeSparsenessMeasure()
        {
            Int32 n = sortedOriginalMSTLengths.Length;
            Int32 n90 = (9 * n) / 10;
            Double sparse = Math.Min(sortedOriginalMSTLengths[n90] / 1000, 1);
            Double t = (Double)totalCount / 500;
            Double correction = .7 + .3 / (1 + t * t);
            return correction * sparse;
        }

        private Double computeStriationMeasure()
        {
            Double numEdges = 0;

            foreach(Edge e in mstEdges)
            {
                Node n1 = e.p1;
                Node n2 = e.p2;
                if (n1.mstDegree == 2 && n2.mstDegree == 2)
                {
                    Edge e1 = getAdjacentMSTEdge(n1, e);
                    Edge e2 = getAdjacentMSTEdge(n2, e);
                    if (cosineOfAdjacentEdges(e, e1, n1) < -.7 && cosineOfAdjacentEdges(e, e2, n2) < -.7)
                        numEdges++;
                }
            }
            return numEdges / (Double)mstEdges.Count;
        }

        private Edge getAdjacentMSTEdge(Node n, Edge e)
        {
            foreach(Edge et in n.neighbors)
            {
                if (et.onMST && e != et)
                {
                    return et;
                }
            }
            return null;
        }

        private Double cosineOfAdjacentEdges(Edge e1, Edge e2, Node n)
        {
            Double v1x = e1.otherNode(n).x - n.x;
            Double v1y = e1.otherNode(n).y - n.y;
            Double v2x = e2.otherNode(n).x - n.x;
            Double v2y = e2.otherNode(n).y - n.y;
            Double v1 = Math.Sqrt(v1x * v1x + v1y * v1y);
            Double v2 = Math.Sqrt(v2x * v2x + v2y * v2y);
            v1x = v1x / v1;
            v1y = v1y / v1;
            v2x = v2x / v2;
            v2y = v2y / v2;
            return v1x * v2x + v1y * v2y;
        }

        private Double computeConvexityMeasure()
        {
            if (hullArea == 0) // poInt32s in general position
                return 1;
            else
            {
                Double t = (Double)totalCount / 500;
                Double correction = .7 + .3 / (1 + t * t);
                Double convexity = alphaArea / hullArea;
                return correction * convexity;
            }
        }

        private Double computeSkinnyMeasure()
        {
            if (alphaPerimeter > 0)
                return 1 - Math.Sqrt(4 * Math.PI * alphaArea) / alphaPerimeter;
            else
                return 1;
        }

        private void computeAlphaArea()
        {
            Double area = 0;       
            foreach(Triangle t in triangles)
            {
                if (t.onComplex)
                {
                    Node p1 = t.anEdge.p1;
                    Node p2 = t.anEdge.p2;
                    Node p3 = t.anEdge.nextE.p2;
                    area += Math.Abs(p1.x * p2.y + p1.y * p3.x + p2.x * p3.y
                            - p3.x * p2.y - p3.y * p1.x - p1.y * p2.x);
                }
            }
            alphaArea = area / 2;
        }

        private void computeHullArea()
        {
            Double area = 0.0;
            foreach(Triangle t in triangles)
            {
                Node p1 = t.anEdge.p1;
                Node p2 = t.anEdge.p2;
                Node p3 = t.anEdge.nextE.p2;
                area += Math.Abs(p1.x * p2.y + p1.y * p3.x + p2.x * p3.y
                        - p3.x * p2.y - p3.y * p1.x - p1.y * p2.x);
            }
            hullArea = area / 2.0;
        }

        private void computeAlphaPerimeter()
        {
            Double sum = 0;
            foreach(Edge e in edges)
            {
                if (e.onShape)
                {
                    sum += e.weight;
                }
            }
            alphaPerimeter = sum;
        }

        private void computeHullPerimeter()
        {
            Double sum = 0;
            Edge e = hullStart;
            do
            {
                sum += e.p1.distToNode(e.p2.x, e.p2.y);
                e = e.nextH;
            } while (!e.isEqual(hullStart));
            hullPerimeter = sum;
        }

        private void setNeighbors()
        {
            foreach(Edge e in edges)
            {
                if (e.isNewEdge(e.p1))
                    e.p1.setNeighbor(e);
                if (e.isNewEdge(e.p2))
                    e.p2.setNeighbor(e);
            }
        }

        private void insert(Int32 px, Int32 py, Int32 count, Int32 id)
        {
            Int32 eid;
            Node nd = new Node(px, py, count, id);
            nodes.Add(nd);
            if (nodes.Count < 3) return;
            if (nodes.Count == 3)    // create the first triangle
            {
                Node p1 = nodes[0];
                Node p2 = nodes[1];
                Node p3 = nodes[2];
                Edge e1 = new Edge(p1, p2);
                if (e1.onSide(p3) == 0)
                {
                    nodes.Remove(nd);
                    return;
                }
                if (e1.onSide(p3) == -1)  // right side
                {
                    p1 = nodes[1];
                    p2 = nodes[0];
                    e1.update(p1, p2);
                }
                Edge e2 = new Edge(p2, p3);
                Edge e3 = new Edge(p3, p1);
                e1.nextH = e2;
                e2.nextH = e3;
                e3.nextH = e1;
                hullStart = e1;
                triangles.Add(new Triangle(edges, e1, e2, e3));
                return;
            }
            actE = edges[0];
            if (actE.onSide(nd) == -1)
            {
                if (actE.invE == null)
                    eid = -1;
                else
                    eid = searchEdge(actE.invE, nd);
            }
            else
                eid = searchEdge(actE, nd);
            if (eid == 0)
            {
                nodes.Remove(nd);
                return;
            }
            if (eid > 0)
                expandTri(actE, nd, eid);   // nd is inside or on a triangle
            else
                expandHull(nd);                // nd is outside convex hull
        }

        private void expandTri(Edge e, Node nd, Int32 type)
        {
            Edge e1 = e;
            Edge e2 = e1.nextE;
            Edge e3 = e2.nextE;
            Node p1 = e1.p1;
            Node p2 = e2.p1;
            Node p3 = e3.p1;
            if (type == 2)
            {   // nd is inside of the triangle
                Edge e10 = new Edge(p1, nd);
                Edge e20 = new Edge(p2, nd);
                Edge e30 = new Edge(p3, nd);
                e.inT.removeEdges(edges);
                triangles.Remove(e.inT);     // remove old triangle
                Edge e100 = e10.makeSymm();
                Edge e200 = e20.makeSymm();
                Edge e300 = e30.makeSymm();
                triangles.Add(new Triangle(edges, e1, e20, e100));
                triangles.Add(new Triangle(edges, e2, e30, e200));
                triangles.Add(new Triangle(edges, e3, e10, e300));
                swapTest(e1);   // swap test for the three new triangles
                swapTest(e2);
                swapTest(e3);
            }
            else
            {   
                // nd is on the edge e
                Edge e4 = e1.invE;
                if (e4 == null || e4.inT == null)
                {          
                    // one triangle involved
                    Edge e30 = new Edge(p3, nd);
                    Edge e02 = new Edge(nd, p2);
                    Edge e10 = new Edge(p1, nd);
                    Edge e03 = e30.makeSymm();
                    //shareEdges(e03,e30);
                    e10.asIndex();
                    e1.mostLeft().nextH = e10;
                    e10.nextH = e02;
                    e02.nextH = e1.nextH;
                    hullStart = e02;
                    triangles.Remove(e1.inT);  // remove oldtriangle and add two new triangles
                    edges.Remove(e1);
                    edges.Add(e10);
                    edges.Add(e02);
                    edges.Add(e30);
                    edges.Add(e03);
                    triangles.Add(new Triangle(e2, e30, e02));
                    triangles.Add(new Triangle(e3, e10, e03));
                    swapTest(e2);   // swap test for the two new triangles
                    swapTest(e3);
                    swapTest(e30);
                }
                else
                {        
                    // two triangle involved
                    Edge e5 = e4.nextE;
                    Edge e6 = e5.nextE;
                    Node p4 = e6.p1;
                    Edge e10 = new Edge(p1, nd);
                    Edge e20 = new Edge(p2, nd);
                    Edge e30 = new Edge(p3, nd);
                    Edge e40 = new Edge(p4, nd);
                    triangles.Remove(e.inT);                   // remove oldtriangle
                    e.inT.removeEdges(edges);
                    triangles.Remove(e4.inT);               // remove old triangle
                    e4.inT.removeEdges(edges);
                    e5.asIndex();   // because e, e4 removed, reset edge sortOrder of node p1 and p2
                    e2.asIndex();
                    triangles.Add(new Triangle(edges, e2, e30, e20.makeSymm()));
                    triangles.Add(new Triangle(edges, e3, e10, e30.makeSymm()));
                    triangles.Add(new Triangle(edges, e5, e40, e10.makeSymm()));
                    triangles.Add(new Triangle(edges, e6, e20, e40.makeSymm()));
                    swapTest(e2);   // swap test for the three new triangles
                    swapTest(e3);
                    swapTest(e5);
                    swapTest(e6);
                    swapTest(e10);
                    swapTest(e20);
                    swapTest(e30);
                    swapTest(e40);
                }
            }
        }

        private void expandHull(Node nd)
        {
            Edge e1, e2, e3 = null, enext;
            Edge e = hullStart;
            Edge comedge = null, lastbe = null;
            while (true)
            {
                enext = e.nextH;
                if (e.onSide(nd) == -1)
                {  
                    // right side
                    if (lastbe != null)
                    {
                        e1 = e.makeSymm();
                        e2 = new Edge(e.p1, nd);
                        e3 = new Edge(nd, e.p2);
                        if (comedge == null)
                        {
                            hullStart = lastbe;
                            lastbe.nextH = e2;
                            lastbe = e2;
                        }
                        else
                            comedge.linkSymm(e2);

                        comedge = e3;
                        triangles.Add(new Triangle(edges, e1, e2, e3));
                        swapTest(e);
                    }
                }
                else
                {
                    if (comedge != null) break;
                    lastbe = e;
                }
                e = enext;
            }

            lastbe.nextH = e3;
            e3.nextH = e;
        }

        private Int32 searchEdge(Edge e, Node nd)
        {
            Int32 f2, f3;
            Edge e0 = null;
            if ((f2 = e.nextE.onSide(nd)) == -1)
            {
                if (e.nextE.invE != null)
                    return searchEdge(e.nextE.invE, nd);
                else
                {
                    actE = e;
                    return -1;
                }
            }
            if (f2 == 0) e0 = e.nextE;
            Edge ee = e.nextE;
            if ((f3 = ee.nextE.onSide(nd)) == -1)
            {
                if (ee.nextE.invE != null)
                    return searchEdge(ee.nextE.invE, nd);
                else
                {
                    actE = ee.nextE;
                    return -1;
                }
            }
            if (f3 == 0) e0 = ee.nextE;
            if (e.onSide(nd) == 0) e0 = e;
            if (e0 != null)
            {
                actE = e0;
                if (e0.nextE.onSide(nd) == 0)
                {
                    actE = e0.nextE;
                    return 0;
                }
                if (e0.nextE.nextE.onSide(nd) == 0) return 0;
                return 1;
            }
            actE = ee;
            return 2;
        }

        private void swapTest(Edge e11)
        {
            Edge e21 = e11.invE;
            if (e21 == null || e21.inT == null) return;
            Edge e12 = e11.nextE;
            Edge e13 = e12.nextE;
            Edge e22 = e21.nextE;
            Edge e23 = e22.nextE;
            if (e11.inT.inCircle(e22.p2) || e21.inT.inCircle(e12.p2))
            {
                e11.update(e22.p2, e12.p2);
                e21.update(e12.p2, e22.p2);
                e11.linkSymm(e21);
                e13.inT.update(e13, e22, e11);
                e23.inT.update(e23, e12, e21);
                e12.asIndex();
                e22.asIndex();
                swapTest(e12);
                swapTest(e22);
                swapTest(e13);
                swapTest(e23);
            }
        }

        private void markHull()
        {
            Edge e = hullStart;
            if (e != null)
                do
                {
                    e.onHull = true;
                    e.p1.onHull = true;
                    e.p2.onHull = true;
                    e = e.nextH;
                } while (!e.isEqual(hullStart));
        }
    }
}
