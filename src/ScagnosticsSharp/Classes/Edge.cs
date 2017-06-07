using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ScagnosticsSharp
{
    internal class Edge
    {
        public Node P1, P2;             // start and end point of the edge
        public Edge InvE = null;        // inverse edge (p2->p1)
        public Edge NextE = null;       // next edge in the triangle in counterclockwise
        public Edge NextH = null;       // convex hull link
        public Triangle InT = null;     // triangle containing this edge
        public Double Weight;
        public Double Coeff_a, Coeff_b, Coeff_c;    // line equation parameters. aX+bY+c=0

        public Boolean OnHull = false;
        public Boolean OnMST = false;
        public Boolean OnShape = false;
        public Boolean OnOutlier = false;

        public Edge(Node p1, Node p2)
        {
            Update(p1, p2);
        }

        public void Update(Node p1, Node p2)
        {
            this.P1 = p1;
            this.P2 = p2;
            Coeff_a = p2.Y - p1.Y;
            Coeff_b = p1.X - p2.X;
            Coeff_c = p2.X * p1.Y - p1.X * p2.Y;
            Weight = Math.Sqrt(Coeff_a * Coeff_a + Coeff_b * Coeff_b);
            AsIndex();
        }

        public Edge MakeSymmEdge()
        {
            Edge e = new Edge(P2, P1);
            LinkSymmEdge(e);
            return e;
        }

        public void LinkSymmEdge(Edge e)
        {
            this.InvE = e;
            if (e != null) e.InvE = this;
        }

        public Int32 OnSide(Node nd)
        {
            Double s = Coeff_a * nd.X + Coeff_b * nd.Y + Coeff_c;
            if (s > 0.0) return 1;
            if (s < 0.0) return -1;
            return 0;
        }

        public void AsIndex()
        {
            P1.AnEdge = this;
        }

        public Edge MostLeft()
        {
            Edge ee, e = this;
            while ((ee = e.NextE.NextE.InvE) != null && ee != this) e = ee;
            return e.NextE.NextE;
        }

        public Edge MostRight()
        {
            Edge ee, e = this;
            while (e.InvE != null && (ee = e.InvE.NextE) != this) e = ee;
            return e;
        }

        public void DeleteSimplex()
        {
            OnShape = false;
            InT.OnComplex = false;
            if (InvE != null)
            {
                InvE.OnShape = false;
                InvE.InT.OnComplex = false;
            }
        }

        public Boolean IsEqual(Edge e)
        {
            return (e.P1.X == this.P1.X) && (e.P2.X == this.P2.X) && (e.P1.Y == this.P1.Y) && (e.P2.Y == this.P2.Y);
        }

        public Boolean IsEquivalent(Edge e)
        {
            return ((e.P1.X == this.P1.X) && (e.P2.X == this.P2.X) && (e.P1.Y == this.P1.Y) && (e.P2.Y == this.P2.Y)) ||
                    ((e.P1.X == this.P2.X) && (e.P1.Y == this.P2.Y) && (e.P2.X == this.P1.X) && (e.P2.Y == this.P1.Y));
        }

        public Node OtherNode(Node n)
        {
            if (n == P1)
                return P2;
            else
                return P1;
        }

        public Boolean IsNewEdge(Node n)
        {
            foreach(Edge e2 in n.Neighbors)
            {
                if (e2.IsEquivalent(this))
                    return false;
            }
            return true;
        }

        public Int32 GetRunts(Double[] maxLength)
        {
            Double cutoff = Weight;
            Double[] maxLength1 = new Double[1];
            Double[] maxLength2 = new Double[1];
            Int32 count1 = P1.GetMSTChildren(cutoff, maxLength1);
            Int32 count2 = P2.GetMSTChildren(cutoff, maxLength2);
            if (count1 < count2)
            {
                maxLength[0] = maxLength1[0];
                return count1;
            }
            else if (count1 == count2)
            {        
                // take more tightly clustered child
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
