using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ScagnosticsSharp
{
    public class Edge
    {
        public Node p1, p2;           // start and end poInt32 of the edge
        protected Edge invE = null;     // inverse edge (p2->p1)
        public Edge nextE = null;    // next edge in the triangle in counterclockwise
        protected Edge nextH = null;    // convex hull link
        public Triangle inT = null;   // triangle containing this edge
        protected Double a, b, c;          // line equation parameters. aX+bY+c=0
        public Double weight;

        protected Boolean onHull = false;
        public Boolean onMST = false;
        protected Boolean onShape = false;
        protected Boolean onOutlier = false;

        protected Edge(Node p1, Node p2)
        {
            update(p1, p2);
        }

        protected void update(Node p1, Node p2)
        {
            this.p1 = p1;
            this.p2 = p2;
            a = p2.y - p1.y;
            b = p1.x - p2.x;
            c = p2.x * p1.y - p1.x * p2.y;
            weight = Math.Sqrt(a * a + b * b);
            asIndex();
        }

        protected Edge makeSymm()
        {
            Edge e = new Edge(p2, p1);
            linkSymm(e);
            return e;
        }

        protected void linkSymm(Edge e)
        {
            this.invE = e;
            if (e != null) e.invE = this;
        }

        protected Int32 onSide(Node nd)
        {
            Double s = a * nd.x + b * nd.y + c;
            if (s > 0.0) return 1;
            if (s < 0.0) return -1;
            return 0;
        }

        protected void asIndex()
        {
            p1.anEdge = this;
        }

        protected Edge mostLeft()
        {
            Edge ee, e = this;
            while ((ee = e.nextE.nextE.invE) != null && ee != this) e = ee;
            return e.nextE.nextE;
        }

        protected Edge mostRight()
        {
            Edge ee, e = this;
            while (e.invE != null && (ee = e.invE.nextE) != this) e = ee;
            return e;
        }

        protected void deleteSimplex()
        {
            onShape = false;
            inT.onComplex = false;
            if (invE != null)
            {
                invE.onShape = false;
                invE.inT.onComplex = false;
            }
        }

        protected Boolean isEqual(Edge e)
        {
            return (e.p1.x == this.p1.x) && (e.p2.x == this.p2.x) && (e.p1.y == this.p1.y) && (e.p2.y == this.p2.y);
        }

        protected Boolean isEquivalent(Edge e)
        {
            return ((e.p1.x == this.p1.x) && (e.p2.x == this.p2.x) && (e.p1.y == this.p1.y) && (e.p2.y == this.p2.y)) ||
                    ((e.p1.x == this.p2.x) && (e.p1.y == this.p2.y) && (e.p2.x == this.p1.x) && (e.p2.y == this.p1.y));
        }

        public Node otherNode(Node n)
        {
            if (n == p1)
                return p2;
            else
                return p1;
        }

        protected Boolean isNewEdge(Node n)
        {
            foreach(Edge e2 in n.neighbors)
            {
                if (e2.isEquivalent(this))
                    return false;
            }
            return true;
        }

        protected Int32 getRunts(Double[] maxLength)
        {
            Double cutoff = weight;
            Double[] maxLength1 = new Double[1];
            Double[] maxLength2 = new Double[1];
            Int32 count1 = p1.getMSTChildren(cutoff, maxLength1);
            Int32 count2 = p2.getMSTChildren(cutoff, maxLength2);
            if (count1 < count2)
            {
                maxLength[0] = maxLength1[0];
                return count1;
            }
            else if (count1 == count2)
            {        // take more tightly clustered child
                if (maxLength1[0] < maxLength2[0])
                    maxLength[0] = maxLength1[0];
                else
                    maxLength[0] = maxLength2[0];
                return count1;
            }
            else
            {
                maxLength[0] = maxLength2[0];
                return count2;
            }
        }
    }
}
