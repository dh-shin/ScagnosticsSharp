using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ScagnosticsSharp
{
    public class Node
    {
        public Int32 x, y;          // coordinate X,Y
        protected Int32 count;        // number of poInt32s aggregated at this node
        public Edge anEdge;     // an edge which starts from this node
        public List<Edge> neighbors;   // nearest Delaunay neighbors list
        protected Boolean onMST;
        protected Boolean onHull = false;
        protected Boolean isVisited = false;
        protected Int32 mstDegree;
        protected Int32 poInt32ID;
        protected Int32 nodeID;

        protected Node(Int32 x, Int32 y, Int32 count, Int32 poInt32ID)
        {
            this.x = x;
            this.y = y;
            this.count = count;
            anEdge = null;
            neighbors = new List<Edge>();
            this.poInt32ID = poInt32ID;
        }

        public Double distToNode(Double px, Double py)
        {
            Double dx = px - x;
            Double dy = py - y;
            return Math.Sqrt(dx * dx + dy * dy);
        }

        protected void setNeighbor(Edge neighbor)
        {
            neighbors.Add(neighbor);
        }

        protected Edge shortestEdge(Boolean mst)
        {
            Edge emin = null;
            if (neighbors != null)
            {
                Double wmin = Double.MaxValue;
                foreach(Edge e in neighbors)
                {
                    if (mst || !e.otherNode(this).onMST)
                    {
                        Double wt = e.weight;
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

        public Int32 getMSTChildren(Double cutoff, Double[] maxLength)
        {
            Int32 count = 0;
            if (isVisited)
                return count;
            isVisited = true;

            foreach(Edge e in neighbors)
            {
                if (e.onMST)
                {
                    if (e.weight < cutoff)
                    {
                        if (!e.otherNode(this).isVisited)
                        {
                            count += e.otherNode(this).getMSTChildren(cutoff, maxLength);
                            Double el = e.weight;
                            if (el > maxLength[0])
                                maxLength[0] = el;
                        }
                    }
                }
            }
            count += this.count; // add count for this node
            return count;
        }
    }
}
