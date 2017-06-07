using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ScagnosticsSharp
{
    internal class Triangle
    {
        public Edge AnEdge;        // an edge of this triangle
        public Double Cx;        // center of circle: X
        public Double Cy;        // center of circle: Y
        public Double Cr;         // radius of circle
        public Boolean OnComplex = true;

        public Triangle(Edge e1, Edge e2, Edge e3)
        {
            Update(e1, e2, e3);
        }

        public Triangle(List<Edge> edges, Edge e1, Edge e2, Edge e3)
        {
            Update(e1, e2, e3);
            edges.Add(e1);
            edges.Add(e2);
            edges.Add(e3);
        }

        public void Update(Edge e1, Edge e2, Edge e3)
        {
            OnComplex = true;
            AnEdge = e1;
            e1.NextE = e2;
            e2.NextE = e3;
            e3.NextE = e1;
            e1.InT = this;
            e2.InT = this;
            e3.InT = this;
            FindCircle();
        }

        public Boolean InCircle(Node nd)
        {
            return nd.DistToNode(Cx, Cy) < Cr;
        }

        public void RemoveEdges(List<Edge> edges)
        {
            edges.Remove(AnEdge);
            edges.Remove(AnEdge.NextE);
            edges.Remove(AnEdge.NextE.NextE);
        }

        private void FindCircle()
        {
            Double x1 = (Double)AnEdge.P1.X;
            Double y1 = (Double)AnEdge.P1.Y;
            Double x2 = (Double)AnEdge.P2.X;
            Double y2 = (Double)AnEdge.P2.Y;
            Double x3 = (Double)AnEdge.NextE.P2.X;
            Double y3 = (Double)AnEdge.NextE.P2.Y;
            Double a = (y2 - y3) * (x2 - x1) - (y2 - y1) * (x2 - x3);
            Double a1 = (x1 + x2) * (x2 - x1) + (y2 - y1) * (y1 + y2);
            Double a2 = (x2 + x3) * (x2 - x3) + (y2 - y3) * (y2 + y3);
            Cx = (a1 * (y2 - y3) - a2 * (y2 - y1)) / a / 2;
            Cy = (a2 * (x2 - x1) - a1 * (x2 - x3)) / a / 2;
            Cr = AnEdge.P1.DistToNode(Cx, Cy);
        }
    }
}
