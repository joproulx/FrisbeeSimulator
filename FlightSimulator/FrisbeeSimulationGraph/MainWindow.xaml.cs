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
using System.Windows.Media.Media3D;
using System.Windows.Navigation;
using System.Windows.Shapes;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics.LinearAlgebra.Generic;
using Petzold.Media3D;

namespace FrisbeeSimulationGraph
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        private double value = 0;
        private Axes2D m_axeXZ;
        private Axes2D m_axeXY;
        private Axes2D m_axeYZ;
        public MainWindow()
        {
            InitializeComponent();

            m_axeXZ = new Axes2D();
            m_axeXY = new Axes2D();
            m_axeYZ = new Axes2D();
            
            panel.Children.Add(m_axeXZ);
            panel.Children.Add(m_axeXY);
            panel.Children.Add(m_axeYZ);
            UpdateLayout();
            
        }

        Matrix<double> GetRotationMatrixX(double angle)
        {
            return new SparseMatrix(new double[,]
                                      {
                                         {1, 0, 0}, 
                                         {0, Math.Cos(angle), -Math.Sin(angle)}, 
                                         {0, Math.Sin(angle), Math.Cos(angle)}
                                     });
        }

        Matrix<double> GetRotationMatrixY  (double angle)
        {
            return new SparseMatrix(new double[,]
                                      {
                                         {Math.Cos(angle), 0, Math.Sin(angle)}, 
                                         {0, 1, 0}, 
                                         { -Math.Sin(angle), 0, Math.Cos(angle)}
                                     });
        }

        Matrix<double> GetRotationMatrixZ(double angle)
        {
            return new SparseMatrix(new double[,]
                                      {
                                         {Math.Cos(angle), -Math.Sin(angle), 0}, 
                                         {Math.Sin(angle), Math.Cos(angle), 0}, 
                                         {0, 0,1}
                                     });
        }


        private void MainWindow_OnKeyDown(object sender, KeyEventArgs e)
        {
        //    Point3D point3D = m_camera.Position;

        //    Vector<double> vector = new SparseVector(3);
        //    vector[0] = point3D.X;
        //    vector[1] = point3D.Y;
        //    vector[2] = point3D.Z;

        //    Vector<double> multiply = null;
        //    if (e.Key == Key.A)
        //    {
        //        multiply = GetRotationMatrixX(Math.PI / 10).Multiply(vector);
        //        multiply[1] = point3D.Y;
        //    }
        //    else if (e.Key == Key.S)
        //    {
        //        vector[0] = 0;
        //        vector[2] = 0;
        //        multiply = GetRotationMatrixY(Math.PI / 10).Multiply(vector);
        //        vector[0] = point3D.X;
        //        vector[2] = point3D.Z;
                
        //    }
        //    else if (e.Key == Key.D)
        //    {
        //        multiply = GetRotationMatrixZ(Math.PI / 10).Multiply(vector);    
        //    }
            
        //    if (multiply != null)
        //    {
        //        m_camera.Position = new Point3D(multiply[0], multiply[1], multiply[2]);
        //        m_camera.LookDirection = new Vector3D(-multiply[0], -multiply[1], -multiply[2]);
        //    }
        }

        private void Window_Loaded(object sender, RoutedEventArgs e)
        {
            
            double theta = 45D*Math.PI/180D, phi = 0D*Math.PI/180D;




            //double x = Math.Sin(angle);
            //double y = -Math.Cos(angle) * Math.Sin(0D);
            //double z = Math.Cos(angle) * Math.Cos(0D);

            //SparseVector vectorSpeed = new SparseVector(new double[]{1D, 1D, 1D });
            //SparseVector vectorC3 = new SparseVector(new double[]{x, y, z });
            //double dotProduct = vectorSpeed.DotProduct(vectorC3);
            //Vector<double> subtract = vectorSpeed.Subtract(vectorC3.Multiply(dotProduct));

            //m_axeXZ.AddPoint(0.5, 0.5);
            //m_axeXZ.AddPoint(subtract[0], subtract[1]);

            Matrix<double> transformation =
                //new SparseMatrix(new double[,]
                //                     {
                //                         {Math.Cos(theta), Math.Sin(theta)*Math.Sin(phi), -Math.Sin(theta)*Math.Cos(phi)}, 
                //                         {0, Math.Cos(phi), Math.Cos(phi)}, 
                //                         {Math.Sin(theta), -Math.Cos(theta)*Math.Cos(phi), Math.Cos(theta)*Math.Cos(phi)}
                //                     });

                new SparseMatrix(new double[,]
                                     {
                                         {0, 0, 0}, 
                                         {0, 0, 0}, 
                                         {Math.Sin(theta), -Math.Cos(theta)*Math.Cos(phi), Math.Cos(theta)*Math.Cos(phi)}
                                     });


            //SparseVector vectorSpeed = new SparseVector(new double[] { 1D, 0D, 0D });

            //Vector<double> multiply = transformation.Multiply(vectorSpeed);

            //DisplayVector(multiply);

            SparseVector vectorC3 = new SparseVector(new double[] { 0D, 0D, 1D });
            SparseVector vectorSpeed = new SparseVector(new double[] { Math.Cos(Math.PI / 4), 0D, Math.Sin(Math.PI / 4) });

            double dotProduct = vectorSpeed.DotProduct(vectorC3);
        }

        private void DisplayVector(Vector<double> multiply)
        {
            m_axeXZ.AddPoint(multiply[0], multiply[2]);
            m_axeXY.AddPoint(multiply[0], multiply[1]);
            m_axeYZ.AddPoint(multiply[1], multiply[2]);
        }
    }
}
