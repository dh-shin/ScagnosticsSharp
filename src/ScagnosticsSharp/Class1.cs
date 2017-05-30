using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ScagnosticsSharp
{
    public class Class1
    {
        private Int32 a;
        private Int32 b;

        public Class1(Int32 a, Int32 b)
        {
            this.a = a;
            this.b = b;
        }

        public void PrintClass1()
        {
            Console.WriteLine(a + b);
        }

        public void PrintClass2()
        {
            Console.WriteLine(a + b + a + b);
        }
    }
}
