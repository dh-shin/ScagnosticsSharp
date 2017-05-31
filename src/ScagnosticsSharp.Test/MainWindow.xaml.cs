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
        private readonly Int32 numBins = 50;    // user setting for number of bins
        private readonly Int32 maxBins = 1000;  // user setting for maximum number of nonempty bins allowed (maxBins >= numBins*numBins)

        private static Int32 VarNum;
        private static Int32 RowNum;
        private static List<List<Double>> rawData;
        private static Double[][] Points;

        public MainWindow()
        {
            InitializeComponent();
            ReadSampleData();

            Scagnostics.LoadJavaRand();
            Double[][] scagnostics = ComputeScagnostics(Points, numBins, maxBins);
            for(Int32 i = 0; i < scagnostics.Length; i++)
            {
                for (Int32 j = 0; j < scagnostics[i].Length ; j++)
                {
                    Console.Write(scagnostics[i][j] + " / ");
                }
                Console.WriteLine("");
            }
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
                    if (isHeader == false)
                    {
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
                        rawData.Add(values);
                    }
                    isHeader = false;
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

        private static Double[][] ComputeScagnostics(Double[][] points, Int32 numBins, Int32 maxBins)
        {
            NormalizePoints(points);
            Int32 nDim = points.Length;
            Int32 numCells = nDim * (nDim - 1) / 2;

            Double[][] scagnostics = new Double[numCells][];
            for(Int32 i = 0; i < numCells; i++)
                scagnostics[i] = new Double[Scagnostics.getNumScagnostics()];

            Int32 k = 0;
            for (Int32 i = 1; i < nDim; i++)
            {
                for (Int32 j = 0; j < i; j++)
                {
                    Scagnostics s = new Scagnostics(points[j], points[i], numBins, maxBins);
                    scagnostics[k] = s.compute();
                    k++;
                }
            }
            return scagnostics;
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
