using System;
using System.Linq;
using System.Collections.Generic;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Media;
using System.Windows.Shapes;

namespace FlightSimulator
{
    class Points
    {
        public List<Point> CartesianPoints { get; set; }
        public Polyline Line { get; set; }

        public Points(Panel parent, Color color)
        {
            CartesianPoints = new List<Point>();
            Line = new Polyline
            {
                Stroke = new SolidColorBrush(color),
                StrokeThickness = 2
            };
            parent.Children.Add(Line);
        }

        public void RefreshLine()
        {

        }
    }

    /// <summary>
    /// Interaction logic for GraphDiagram.xaml
    /// </summary>
    public partial class GraphDiagram : UserControl
    {

        private readonly Dictionary<Color, Points> m_points;
        private Polyline m_lineX;
        private Polyline m_lineY;

        public GraphDiagram()
        {
            InitializeComponent();

            MaxValueX = double.MinValue;
            MaxValueY = double.MinValue;
            MinValueX = double.MaxValue;
            MinValueY = double.MaxValue;
            m_points = new Dictionary<Color, Points>();

            m_lineX = new Polyline
            {
                Stroke = new SolidColorBrush(Colors.Black),
                StrokeThickness = 2
            };
            m_lineY = new Polyline
            {
                Stroke = new SolidColorBrush(Colors.Black),
                StrokeThickness = 2
            };

            canvas.Children.Add(m_lineX);
            canvas.Children.Add(m_lineY);
        }

        protected double MaxValueX { get; set; }

        protected double MaxValueY { get; set; }

        protected double MinValueX { get; set; }

        protected double MinValueY { get; set; }


        public void AddPoint(double x, double y, Color color)
        {
            Points pointList = GetPoints(color);

            pointList.CartesianPoints.Add(new Point(x, y));

            if (x > MaxValueX)
            {
                MaxValueX = x;
            }

            if (y > MaxValueY)
            {
                MaxValueY = y;
            }

            if (y < MinValueY)
            {
                MinValueY = y;
            }

            if (x < MinValueX)
            {
                MinValueX = x;
            }
        }

        private Points GetPoints(Color color)
        {
            Points points;
            if (!m_points.TryGetValue(color, out points))
            {
                points = new Points(canvas, color);
                m_points.Add(color, points);
            }

            return points;
        }

        public void Reset()
        {
            foreach (var points in m_points.Values)
            {
                points.CartesianPoints.Clear();
            }
        }

        public void Refresh()
        {
            foreach (Points points in m_points.Values)
            {
                points.Line.Points.Clear();
                Points points1 = points;
                points.CartesianPoints.ForEach(p => points1.Line.Points.Add(Adjust(p)));
            }

            m_lineX.Points.Clear();
            m_lineY.Points.Clear();

            m_lineX.Points.Add(Adjust(new Point(0, MinValueY)));
            m_lineX.Points.Add(Adjust(new Point(0, MaxValueY)));
            m_lineX.Points.Add(Adjust(new Point(MinValueX, 0)));
            m_lineX.Points.Add(Adjust(new Point(MaxValueX, 0)));
        }

        private Point Adjust(Point point)
        {
            var newPoint = new Point
                               {
                                   X = point.X - MinValueX,
                                   Y = point.Y - MinValueY
                               };

            double width = MaxValueX - MinValueX;
            double height = MaxValueY - MinValueY;

            double graphWidth = canvas.ActualWidth;
            double graphHeight = canvas.ActualHeight;

            newPoint.X = (newPoint.X * graphWidth) / width;
            newPoint.Y = graphHeight - ((newPoint.Y * graphHeight) / height);

            return newPoint;
        }

        //public void AddPoints(List<Point> simulate)
        //{
        //    simulate.ForEach(p => AddPoint(p.X, p.Y));
        //}
    }
}