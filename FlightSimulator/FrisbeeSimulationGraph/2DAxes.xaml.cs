using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics.LinearAlgebra.Generic;

namespace FrisbeeSimulationGraph
{
    /// <summary>
    /// Interaction logic for _2DAxes.xaml
    /// </summary>
    public partial class Axes2D : UserControl
    {
        public double MaxX { get; set; }
        public double MaxY { get; set; }

        private readonly List<Point> m_points;

        public Axes2D()
        {
            m_points = new List<Point>();

            MaxX = double.MinValue;
            MaxY = double.MinValue;
            InitializeComponent();

        }


        public void AddPoint(double x, double y)
        {
            m_points.Add(new Point(x, y));
            if (Math.Abs(x) > MaxX)
            {
                MaxX = Math.Abs(x);
            }

            if (Math.Abs(y) > MaxY)
            {
                MaxY = Math.Abs(y);
            }
        }

        public void DisplayPoints()
        {
            foreach (Point point in m_points)
            {
                double newX = 0.5D + point.X / (2 * MaxX);
                double newY = 0.5D - point.Y / (2 * MaxY);

                Line line = new Line
                {
                    X1 = 0.5d*ActualWidth,
                    Y1 = 0.5d*ActualHeight,
                    X2 = newX*ActualWidth,
                    Y2 = newY*ActualHeight,
                    Stroke = new SolidColorBrush(Colors.Red)
                };

                m_layout.Children.Add(line);
            }
        }

    }
}
