using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ScagnosticsSharp
{
    internal class Node
    {
        public Int32 X, Y;              // coordinate X, Y
        public Int32 Count;             // number of points aggregated at this node
        public Edge AnEdge;             // an edge which starts from this node
        public List<Edge> Neighbors;    // nearest Delaunay neighbors list
        public Boolean OnMST;
        public Boolean OnHull = false;
        public Boolean IsVisited = false;
        public Int32 MstDegree;
        public Int32 PointID;

        public Node(Int32 x, Int32 y, Int32 count, Int32 pointID)
        {
            X = x;
            Y = y;
            Count = count;
            AnEdge = null;
            Neighbors = new List<Edge>();
            PointID = pointID;
        }

        public Double DistToNode(Double px, Double py)
        {
            Double dx = px - X;
            Double dy = py - Y;
            return Math.Sqrt(dx * dx + dy * dy);
        }

        public void SetNeighbor(Edge neighbor)
        {
            Neighbors.Add(neighbor);
        }

        public Edge ShortestEdge(Boolean mst)
        {
            Edge emin = null;
            if (Neighbors != null)
            {
                Double wmin = Double.MaxValue;
                foreach(Edge e in Neighbors)
                {
                    if (mst || !e.OtherNode(this).OnMST)
                    {
                        Double wt = e.Weight;
                        if (wt < wmin)
                        {
                            wmin = wt;
                            emin = e;
                        }
                    }
                }
            }
            return emin;
        }

        public Int32 GetMSTChildren(Double cutoff, Double[] maxLength)
        {
            Int32 count = 0;
            if (IsVisited)
                return count;
            IsVisited = true;

            foreach(Edge e in Neighbors)
            {
                if (e.OnMST)
                {
                    if (e.Weight < cutoff)
                    {
                        if (!e.OtherNode(this).IsVisited)
                        {
                            count += e.OtherNode(this).GetMSTChildren(cutoff, maxLength);
                            Double el = e.Weight;
                            if (el > maxLength[0])
                                maxLength[0] = el;
                        }
                    }
                }
            }
            count += Count; // add count for this node
            return count;
        }
    }
}
