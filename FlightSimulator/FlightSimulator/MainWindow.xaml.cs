using System;
using System.Collections.Generic;
using System.Windows;
using System.Windows.Media;

namespace FlightSimulator
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        double[] m_y0Initials = new double[]
                                  {
                                      -9.03E-01,
                                      -6.33E-01,
                                      -9.13E-01,
                                      1.34E+01,
                                      -4.11E-01,
                                      1.12E-03,
                                      -7.11E-02,
                                      2.11E-01,
                                      -1.49E+01,
                                      -1.48E+00,
                                      5.43E+01,
                                      5.03E+00
                                  };

        private double[] m_y0;


        private readonly Frisbee m_frisbee = new Frisbee();

        public MainWindow()
        {
            InitializeComponent();

            m_y0 = new double[m_y0Initials.Length];

            Reset();


        }

        private void Reset()
        {
            for (int i = 0; i < m_y0Initials.Length; i++)
            {
                m_y0[i] = m_y0Initials[i];
            }
        }

        private void OnGraphLoaded(object sender, RoutedEventArgs e)
        {


            ExecuteSimulation(m_y0);

            /* List<Point> simulate = m_frisbee.Simulate(1, 14, 0, 90, 0.001);
            m_graph.AddPoints(simulate);
            m_graph.Refresh();*/
        }

        private void ExecuteSimulation(double[] y0)
        {
            m_graph.Reset();
            try
            {
                
                

                //double tfinal = 1.46; //% length of flight
                double tfinal = 5; //% length of flight
                double nsteps = 292;// % number of time steps for data

                double span = tfinal/nsteps;

                double[] x = new double[(int)nsteps];
                for (int i = 0; i < nsteps; i++)
                {
                    x[i] = span*i;
                }

                double[,] y = m_frisbee.Ode(y0, x); 


                //using (System.IO.StreamWriter file = new System.IO.StreamWriter(@"D:\temp\test.txt"))
                {
                    for (int i = 0; i < y.GetLength(0); i++)
                    {
                        for (int j = 0; j < y.GetLength(1); j++)
                        {
                            if (j == 1)
                            {
                                m_graph.AddPoint(y[i, 0], y[i, j], Colors.Blue);
                            }

                            if (j == 2)
                            {
                                m_graph.AddPoint(y[i, 0], y[i, j], Colors.Red);
                            }

                            if (j == 3)
                            {
                                m_graph.AddPoint(y[i, 0], y[i, j], Colors.Green);
                            }

                            //file.Write(y[i, j] + "\t");
                        }
                        //file.WriteLine("");
                    }

                    
                }

                Dispatcher.BeginInvoke(new Action(() => m_graph.Refresh()));


                //m_frisbee.Simulate(new Frisbee.SimulationState
                //{
                //    VX = 1.34E+01,
                //    VY = -4.11E-01,
                //    VZ = 1.12E-03,
                //    Phi = -7.11E-02,
                //    Theta = 2.11E-01,
                //    PhiDot = -1.49E+01,
                //    ThetaDot = -1.48E+00,
                //    GammaDot = 5.43E+01,
                //});
            }
            catch (Exception ex)
            {

                MessageBox.Show(ex.ToString());
            }
        }

        private void OnGraphSizeChanged(object sender, SizeChangedEventArgs e)
        {
            m_graph.Refresh();
        }

        private void buttonRefresh_Click(object sender, RoutedEventArgs e)
        {
            m_y0[3] = sliderVx.Value;
            m_y0[4] = -sliderVy.Value;
            m_y0[5] = sliderVz.Value;
            m_y0[6] = sliderPhi.Value;
            m_y0[7] = sliderTheta.Value;
            m_y0[11] = sliderGamma.Value;
            ExecuteSimulation(m_y0);
            m_graph.Refresh();
        }   

        private void buttonReset_Click(object sender, RoutedEventArgs e)
        {
            Reset();
        }
    }
}