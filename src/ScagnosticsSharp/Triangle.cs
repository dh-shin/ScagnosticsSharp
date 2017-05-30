using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ScagnosticsSharp
{
    public class Triangle
    {
        public Edge anEdge;        // an edge of this triangle
        protected Double c_cx;        // center of circle: X
        protected Double c_cy;        // center of circle: Y
        protected Double c_r;         // radius of circle
        public Boolean onComplex = true;

        public Triangle(Edge e1, Edge e2, Edge e3)
        {
            update(e1, e2, e3);
        }

        public Triangle(List<Edge> edges, Edge e1, Edge e2, Edge e3)
        {
            update(e1, e2, e3);
            edges.Add(e1);
            edges.Add(e2);
            edges.Add(e3);
        }

        public void update(Edge e1, Edge e2, Edge e3)
        {
            onComplex = true;
            anEdge = e1;
            e1.nextE = e2;
            e2.nextE = e3;
            e3.nextE = e1;
            e1.inT = this;
            e2.inT = this;
            e3.inT = this;
            findCircle();
        }

        public Boolean inCircle(Node nd)
        {
            return nd.distToNode(c_cx, c_cy) < c_r;
        }

        public void removeEdges(List<Edge> edges)
        {
            edges.Remove(anEdge);
            edges.Remove(anEdge.nextE);
            edges.Remove(anEdge.nextE.nextE);
        }

        protected void findCircle()
        {
            Double x1 = (Double)anEdge.p1.x;
            Double y1 = (Double)anEdge.p1.y;
            Double x2 = (Double)anEdge.p2.x;
            Double y2 = (Double)anEdge.p2.y;
            Double x3 = (Double)anEdge.nextE.p2.x;
            Double y3 = (Double)anEdge.nextE.p2.y;
            Double a = (y2 - y3) * (x2 - x1) - (y2 - y1) * (x2 - x3);
            Double a1 = (x1 + x2) * (x2 - x1) + (y2 - y1) * (y1 + y2);
            Double a2 = (x2 + x3) * (x2 - x3) + (y2 - y3) * (y2 + y3);
            c_cx = (a1 * (y2 - y3) - a2 * (y2 - y1)) / a / 2;
            c_cy = (a2 * (x2 - x1) - a1 * (x2 - x3)) / a / 2;
            c_r = anEdge.p1.distToNode(c_cx, c_cy);
        }
    }
}
