using System;
using System.Collections.Generic;
using System.Linq;
using System.IO;
using System.Reflection;
using System.Windows;

namespace ScagnosticsSharp.Test
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        private static Int32 VarNum;
        private static Int32 RowNum;
        private static List<List<Double>> rawData;
        private static String[] Headers;
        private static Double[][] Points;

        public MainWindow()
        {
            InitializeComponent();

            ReadSampleData();

            Scagnostics.LoadJavaRandomNumber();

            Console.WriteLine("================= RUN SCAGNOSTICS =================");

            NormalizePoints(Points);

            Int32 nDim = Points.Length;
            Int32 numCells = nDim * (nDim - 1) / 2;

            String[][] vars = new String[numCells][];
            for (Int32 i = 0; i < numCells; i++)
                vars[i] = new String[2];

            Double[][] results = new Double[numCells][];
            for (Int32 i = 0; i < numCells; i++)
                results[i] = new Double[Scagnostics.GetNumScagnostics()];

            Int32 k = 0;
            for (Int32 i = 1; i < nDim; i++)
            {
                for (Int32 j = 0; j < i; j++)
                {
                    Scagnostics scagnostics = new Scagnostics(Points[j], Points[i]);
                    results[k] = scagnostics.Compute();
                    
                    vars[k][0] = Headers[j];
                    vars[k][1] = Headers[i];

                    k++;
                }
            }

            String[] measures = Scagnostics.GetScagnosticsLabels();
            for (Int32 i = 0; i < measures.Length; i++)
                Console.Write(measures[i] + " / ");
            Console.WriteLine("");

            for (Int32 i = 0; i < results.Length; i++)
            {
                Console.Write(vars[i][0] + " - " + vars[i][1] + " : ");

                for (Int32 j = 0; j < results[i].Length; j++)
                    Console.Write(results[i][j] + " / ");
                Console.WriteLine("");
            }

            Console.WriteLine("===================================================");
        }

        private static void ReadSampleData()
        {   
            var assembly = Assembly.GetExecutingAssembly();
            var resourceName = "ScagnosticsSharp.Test.Resources.baseball.csv";
            using (Stream stream = assembly.GetManifestResourceStream(resourceName))
            using (StreamReader reader = new StreamReader(stream))
            {
                rawData = new List<List<Double>>();
                Boolean isHeader = true;
                Char[] delimChars = { ',' };
                while (reader.EndOfStream == false)
                {
                    String result = reader.ReadLine();
                    String[] tokenized = result.Split(delimChars);

                    if (isHeader == true)
                    {
                        Headers = tokenized;
                        isHeader = false;
                    }
                    else
                    {
                        List<Double> values = new List<Double>();
                        for (Int32 i = 0; i < tokenized.Length; i++)
                        {
                            Double val;
                            if (Double.TryParse(tokenized[i], out val) == true)
                            {
                                values.Add(val);
                            }
                        }
                        rawData.Add(values);
                    }
                }
            }

            VarNum = rawData.First().Count;
            RowNum = rawData.Count;

            Points = new Double[VarNum][];
            for (Int32 i = 0; i < VarNum; i++)
                Points[i] = new Double[RowNum];

            for (Int32 i = 0; i < VarNum; i++)
            {
                for (Int32 j = 0; j < RowNum; j++)
                {
                    Points[i][j] = rawData[j][i];
                }
            }
        }

        private static Double[][] ComputeScagnostics(Double[][] points)
        {
            NormalizePoints(points);
            Int32 nDim = points.Length;
            Int32 numCells = nDim * (nDim - 1) / 2;

            Double[][] results = new Double[numCells][];
            for(Int32 i = 0; i < numCells; i++)
                results[i] = new Double[Scagnostics.GetNumScagnostics()];

            Int32 k = 0;
            for (Int32 i = 1; i < nDim; i++)
            {
                for (Int32 j = 0; j < i; j++)
                {
                    Scagnostics scagnostics = new Scagnostics(points[j], points[i]);
                    results[k] = scagnostics.Compute();
                    k++;
                }
            }
            return results;
        }

        private static void NormalizePoints(Double[][] points)
        {
            Double[] min = new Double[points.Length];
            Double[] max = new Double[points.Length];
            for (Int32 i = 0; i < points.Length; i++)
                for (Int32 j = 0; j < points[0].Length; j++)
                {
                    if (j == 0)
                        max[i] = min[i] = points[i][0];
                    else if (min[i] > points[i][j])
                        min[i] = points[i][j];
                    else if (max[i] < points[i][j])
                        max[i] = points[i][j];
                }
            for (Int32 i = 0; i < points.Length; i++)
                for (Int32 j = 0; j < points[0].Length; j++)
                    points[i][j] = (points[i][j] - min[i]) / (max[i] - min[i]);
        }
    }
}
