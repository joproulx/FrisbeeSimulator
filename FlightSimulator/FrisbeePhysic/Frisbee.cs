using System;
using System.Collections.Generic;
using System.Drawing;
using System.Numerics;
using System.Windows.Media.Media3D;
using DotNumerics.ODE;
using MathNet.Numerics.Interpolation.Algorithms;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics.LinearAlgebra.Double.Factorization;
using MathNet.Numerics.LinearAlgebra.Generic;

namespace FlightSimulator
{
    public class Frisbee
    {
        static readonly List<double> m_xCL = new List<double>
                                    {
                                        -0.1745,
                                        -0.05236,
                                        0,
                                        0.08727,
                                        0.17453,
                                        0.26180,
                                        0.34907,
                                        0.43633,
                                        0.52360,
                                    };

        static readonly List<double> m_yCL = new List<double>
                                    {
                                        -0.2250, 
                                        0, 
                                        0.150, 
                                        0.4500, 
                                        0.7250, 
                                        0.9750, 
                                        1.2000, 
                                        1.4500,
                                        1.6750, 
                                    };


        static readonly List<double> m_xCD = new List<double>
                                    {
                                        -0.1745,   
                                        -0.05236,  
                                        0,         
                                        0.08727,   
                                        0.1745,    
                                        0.26180,   
                                        0.3491,    
                                        0.4363,    
                                        0.5236    
                                    };

        static readonly List<double> m_yCD = new List<double>
                                    {
                                        0.1500,
                                        0.0800,
                                        0.1000,
                                        0.1500,
                                        0.2600,
                                        0.3900,
                                        0.5700,
                                        0.7500,
                                        0.9200,
                                    };
        
        static readonly List<double> m_xCM = new List<double>
                                    {
                                        -0.174532925,
                                        -0.087266463,
                                        -0.052359878,
                                        0,           
                                        0.052359878, 
                                        0.104719755, 
                                        0.157079633, 
                                        0.20943951,  
                                    };
        
         static readonly List<double> m_yCM = new List<double>
                                    {
                                        -0.0380,
                                        -0.0220,
                                        -0.0140,
                                        -0.0060,
                                        -0.0060,
                                        -0.0020,
                                        0.0000, 
                                        0.0100, 
                                    };



        private double x;
        //The x position of the frisbee.
        private double y;
        //The y position of the frisbee.
        private double vx;
        //The x velocity of the frisbee.
        private double vy;
        //The y velocity of the frisbee.
        private const double g = -9.81;

        //The acceleration of gravity (m/s^2).
        private const double m = 0.175;
        //The mass of a standard frisbee in kilograms.
        private const double Rho = 1.23;
        //The density of air in kg/m^3.
        private const double Area = 0.0568;
        //The area of a standard frisbee.
        private const double CL0 = 0.1;
        //The lift coefficient at alpha = 0.
        private const double CLA = 1.4;
        //The lift coefficient dependent on alpha.
        private const double CD0 = 0.08;
        //The drag coefficent at alpha = 0.
        private const double CDA = 2.72;
        //The drag coefficient dependent on alpha.
        private const double ALPHA0 = -4;

        private const double Ia = 0.002352; //% moment of inertia about the spinning axis
        private const double Id = 0.001219; //% moment of inertia about the planar axis'


        /**
        * A method that uses Euler’s method to simulate the flight of a frisbee in
        * two dimensions, distance and height (x and y, respectively).
        *
        */
        //public List<Point> Simulate(double y0, double vx0, double vy0, double alpha, double deltaT)
        //{
        //    List<Point> points = new List<Point>();
        //    //Calculation of the lift coefficient using the relationship given
        //    //by S. A. Hummel.
        //    double cl = CL0 + CLA * alpha * Math.PI / 180;
        //    //Calculation of the drag coefficient (for Prantl’s relationship)
        //    //using the relationship given by S. A. Hummel.
        //    double cd = CD0 + CDA * Math.Pow((alpha - ALPHA0) * Math.PI / 180, 2);
        //    //Initial position x = 0.
        //    x = 0;
        //    //Initial position y = y0.
        //    y = y0;
        //    //Initial x velocity vx = vx0.
        //    vx = vx0;
        //    //Initial y velocity vy = vy0.
        //    vy = vy0;
        //    try
        //    {

        //        //A PrintWriter object to write the output to a spreadsheet.
        //        //A loop index to monitor the simulation steps.
        //        int k = 0;
        //        //A while loop that performs iterations until the y position
        //        //reaches zero (i.e. the frisbee hits the ground).
        //        while (y > 0)
        //        {
        //            //The change in velocity in the y direction obtained setting the
        //            //net force equal to the sum of the gravitational force and the
        //            //lift force and solving for delta v.
        //            double deltavy = (Rho * Math.Pow(vx, 2) * Area * cl / 2 / m + g) * deltaT;
        //            //The change in velocity in the x direction, obtained by
        //            //solving the force equation for delta v. (The only force
        //            //present is the drag force).
        //            double deltavx = -Rho * Math.Pow(vx, 2) * Area * cd * deltaT;
        //            //The new positions and velocities are calculated using
        //            //simple introductory mechanics.
        //            vx = vx + deltavx;
        //            vy = vy + deltavy;
        //            x = x + vx * deltaT;
        //            y = y + vy * deltaT;
        //            //Only the output from every tenth iteration will be sent
        //            //to the spreadsheet so as to decrease the number of data points.
        //            if (k / 10 == 0)
        //            {
        //                points.Add(new Point(x, y));
        //            }
        //            k++;
        //        }
        //    }
        //    catch (Exception)
        //    {
        //        Console.WriteLine(@"Error, file frisbee.csv is in use.");
        //    }
        //    return points;
        //}

        public List<Point> Simulate(double y0, double vx0, double vy0, double alpha, double deltaT)
        {
            List<Point> points = new List<Point>();
            //Calculation of the lift coefficient using the relationship given
            //by S. A. Hummel.
            double cl = CL0 + CLA * alpha * Math.PI / 180;
            //Calculation of the drag coefficient (for Prantl’s relationship)
            //using the relationship given by S. A. Hummel.
            double cd = CD0 + CDA * Math.Pow((alpha - ALPHA0) * Math.PI / 180, 2);
            //Initial position x = 0.
            x = 0;
            //Initial position y = y0.
            y = y0;
            //Initial x velocity vx = vx0.
            vx = vx0;
            //Initial y velocity vy = vy0.
            vy = vy0;
            try
            {
                

      
            }
            catch (Exception)
            {
                Console.WriteLine(@"Error, file frisbee.csv is in use.");
            }
            return points;
        }


        public class SimulationState
        {
            public double VX { get; set; }
            public double VY { get; set; }
            public double VZ { get; set; }
            public double Phi { get; set; }
            public double Theta { get; set; }
            public double Gamma { get; set; }
            public double SinTheta { get { return Math.Sin(Theta); } }
            public double CosTheta { get { return Math.Cos(Theta); } }
            public double SinPhi { get { return Math.Sin(Phi); } }
            public double CosPhi { get { return Math.Cos(Phi); } }
            
            public double PhiDot { get; set; }
            public double ThetaDot { get; set; }
            public double GammaDot { get; set; }
        }

        public double[,] Ode(double[] y0, double[] x)
        {
            OdeExplicitRungeKutta45 RK45 = new OdeExplicitRungeKutta45(Simulate, 12)
            {
                RelTol = 0.000001,
                AbsTol = 0.000001
            };

            return RK45.Solve(y0, x);
        }
        
        public double[] Simulate(double t, double[] values)
        {
            SimulationState st = new SimulationState
                                     {
                                         VX = values[3],
                                         VY = values[4],
                                         VZ = values[5],
                                         Phi = values[6],
                                         Theta = values[7],
                                         PhiDot = values[8],
                                         ThetaDot = values[9],
                                         GammaDot = values[10],
                                         Gamma = values[11],
                                     };
            
            //double CLo, CLa, CDo, CDa, CMo, CMa;
            //double CL_data, CD_data, CM_data;
            
            double[] CRr_rad = new[] { -0.0873, -0.0698, -0.0524, -0.0349, -0.0175, 0.0000, 0.0175, 0.0349, 0.0524, 0.0698, 0.0873, 0.1047, 0.1222, 0.1396, 0.1571, 0.1745, 0.1920, 0.2094, 0.2269, 0.2443, 0.2618, 0.5236 };

            double[] CRr_AdvR = new[] {2,1.04,0.69,0.35,0.17,0};

            double[,] CRr_data = new [,]
                                     {
                                         { -0.0172, -0.0192, -0.018, -0.0192, -0.018, -0.0172, -0.0172, -0.0168, -0.0188, -0.0164, -0.0136, -0.01, -0.0104, -0.0108, -0.0084, -0.008, -0.008, -0.006, -0.0048, -0.0064, -0.008, -0.003 },
                                         { -0.0112, -0.0132, -0.012, -0.0132, -0.012, -0.0112, -0.0112, -0.0108, -0.0128, -0.0104, -0.0096, -0.0068, -0.0072, -0.0076, -0.0052, -0.0048, -0.0048, -0.0028, -0.0032, -0.0048, -0.0064, -0.003 },
                                         { -0.0056, -0.0064, -0.0064, -0.0068, -0.0064, -0.0064, -0.0052, -0.0064, -0.0028, -0.0028, -0.004, -0.002, -0.004, -0.002, -0.0016, 0, 0, 0, 0, -0.002, -0.0048, -0.003 },
                                         { -0.0012, -0.0016, -0.0004, -0.0028, -0.0016, -0.0016, -0.0004, 0.0004, 0.0004, 0.0008, 0.0004, 0.0008, 0.0012, 0.0008, 0.002, 0.0028, 0.0032, 0.0024, 0.0028, 0.0004, -0.0012, -0.003 },
                                         { -0.0012, -0.0012, -0.0016, -0.0016, -0.0012, -0.0004, 0.0004, 0.0008, 0.0008, 0.0016, 0.0004, 0.002, 0.0004, 0.0016, 0.002, 0.002, 0.002, 0.0012, 0.0012, 0, -0.0012, -0.003 },
                                         { -0.0012, -0.0012, -0.0004, -0.0008, -0.0008, -0.0008, 0.0004, 0.0004, 0.0004, 0.0008, 0.0004, 0.0008, -0.0004, 0, 0, 0.0004, 0, 0, 0.0004, -0.002, -0.0012, -0.003 }
                                     };
            
            const double CMq = -1.44E-02;
            const double CRp = -1.25E-02;
            const double CNr = -3.41E-05;
            
            double diameter = 2*Math.Sqrt(Area/Math.PI);
            
            //  Rotation matrix: http://s-mat-pcs.oulu.fi/~mpa/matreng/eem1_3-7.htm
            //                           y      
            //  ------------------> x                        ^ 
            //  |\                       |    
            //  | \                      |   
            //  |  \                     |    
            //  |   \ theta = pitch      |   gamma = yaw    
            //  |                        |         
            //  v                        --------------------> x        
            //  z
            //  z
            //  ^
            //  |
            //  |
            //  |
            //  |  phi = roll
            //  |
            //  ------------------> y
            //
            // 3D homogenous transformation matrix
            // 
            // g = gamma = yaw
            // t = theta = pitch
            // p = phi   = roll
            // 
            // http://en.wikipedia.org/wiki/Rotation_matrix
            // http://www.gregslabaugh.name/publications/euler.pdf
            //       --                                                                                   --              
            //      | cos(g)*cos(t), cos(g)*sin(t)*sin(p)-sin(g)*cos(p), cos(g)*sin(t)*cos(p)+sin(g)*sin(p) | 
            //      |                                                                                       |
            //  T = | sin(g)*cos(t), sin(g)*sin(t)*sin(p)-cos(g)*cos(p), sin(g)*sin(t)*cos(p)+cos(g)*sin(p) |
            //      |                                                                                       |
            //      |    -sin(t)   ,          cos(t)*sin(p)            ,          cos(t)*cos(p)             |                
            //       --                                                                                   --
            //
            // With g = yaw = 0 and sin(t) = -sin(t) since z is positive downward 
            //
            //       --                                      --              
            //      | cos(t)  , sin(t)*sin(p)  , sin(t)*cos(p) | 
            //      |                                          |
            //  T = |   0     ,     cos(p)     ,     sin(p)    |
            //      |                                          |
            //      | -sin(t) , cos(t)*sin(p)  , cos(t)*cos(p) |                
            //       --                                      --
            
            
            
            Matrix<double> transformation = new SparseMatrix(new [,]
            {
                {st.CosTheta, st.SinTheta*st.SinPhi, -st.SinTheta*st.CosPhi}, 
                {0, st.CosPhi, st.SinPhi}, 
                {st.SinTheta, -st.CosTheta*st.SinPhi, st.CosTheta*st.CosPhi}
            });

            // Eigenvector & eigenvalue
            //       --
            //      | x1
            //  X = | x2
            //      | x3
            //       --
            //
            //       --
            //      | a11, a12, a13
            //  A = | a21, a22, a23
            //      | a31, a32, a33
            //       --
            //
            //
            // Usually, the multiplication of a vector x by a square matrix A changes both the magnitude and the direction 
            // of the vector upon which it acts; but in the special case where it changes only the scale (magnitude) of the 
            // vector and leaves the direction unchanged, or switches the vector to the opposite direction, then that vector 
            // is called an eigenvector of that matrix (the term "eigenvector" is meaningless except in relation to some 
            // particular matrix). When multiplied by a matrix, each eigenvector of that matrix changes its magnitude by a 
            // factor, called the eigenvalue corresponding to that eigenvector. 
            //
            //
            // See local frame vs global frame: 
            //
            
            Evd evd = new UserEvd(transformation);
            //Matrix<double> eigenVectors = evd.EigenVectors();
            Vector<Complex> temp = evd.EigenValues();
            Vector<double> eigenValues = new SparseVector(3);

            eigenValues[0] = temp[0].Real;
            eigenValues[1] = temp[1].Real;
            eigenValues[1] = temp[1].Real;

            //eigenValues.Norm

            //
            // If you have Theta and Phi = 0 you have a transformation matrix like this:
            // 
            // | 1, 0, 0 |
            // | 0, 1, 0 |
            // | 0, 0, 1 |
            // 
            // So each row represents the rotated X, Y or Z axis expressed as N-Frame coordinates. In this case,
            // there is no rotation so you have the axis (1,0,0), (0,1,0), (0,0,1). 
            // For example, the first row represents the X Axis after the rotation. Since the rotation is 0,
            // the X axis is a vector (1,0,0) in the N-Frame.
            //
            //
            //
            //
            //SparseVector c1 = new SparseVector(transformation.Row(0));
            //SparseVector c2 = new SparseVector(transformation.Row(1));
            SparseVector c3 = new SparseVector(transformation.Row(2));

            SparseVector velocity = new SparseVector(new []{st.VX, st.VY, st.VZ});
            double velocityMagnitude = velocity.Norm(2);

            double velocityC3 = velocity.DotProduct(c3);

            Vector<double> vp = velocity.Subtract(c3.Multiply(velocityC3));
            double vpMagnitude = vp.Norm(2);

            
            double alpha = Math.Atan(velocityC3/vp.Norm(2));
            double adp = Area*Rho*velocityMagnitude*velocityMagnitude/2;

            Vector<double> unitVelocity = velocity.Divide(velocityMagnitude);
            Vector<double> unitVp = vp.Divide(vpMagnitude);

            //c3.
            Vector<double> unitLat = ConvertVector(Vector3D.CrossProduct(ConvertVector(c3), ConvertVector(unitVp)));

            Matrix<double> omegaD_N_inC = new SparseMatrix(new [,]{{st.PhiDot*st.CosTheta, st.ThetaDot, st.PhiDot*st.SinTheta+st.GammaDot}}); // expressed in c1,c2,c3
            Vector<double> omegaD_N_inN = transformation.Transpose().Multiply(omegaD_N_inC.Transpose()).Column(0); // expressed in c1,c2,c3
            double omegaVp = omegaD_N_inN.DotProduct(unitVp);
            double omegaLat = omegaD_N_inN.DotProduct(unitLat);
            double omegaSpin = omegaD_N_inN.DotProduct(c3);

            double aDvR = diameter*omegaSpin/2/vpMagnitude;



            LinearSplineInterpolation interpolation = new LinearSplineInterpolation(m_xCL, m_yCL);
            double CL = interpolation.Interpolate(alpha);

            interpolation = new LinearSplineInterpolation(m_xCD, m_yCD);
            double CD = interpolation.Interpolate(alpha);

            interpolation = new LinearSplineInterpolation(m_xCM, m_yCM);
            double CM = interpolation.Interpolate(alpha);


            alglib.spline2d.spline2dinterpolant interpolant = new alglib.spline2d.spline2dinterpolant();
            alglib.spline2d.spline2dbuildbilinear(CRr_rad, CRr_AdvR, CRr_data, 6, 22, interpolant);
            double CRr = alglib.spline2d.spline2dcalc(interpolant, alpha, aDvR);

            Vector<double> mvp = unitVp.Multiply(adp * diameter * (CRr * diameter * omegaSpin / 2 / velocityMagnitude + CRp * omegaVp)); // Roll moment, expressed in N
            
            double lift = CL*adp;
            double drag = CD*adp;
            
            Vector<double> unitLift = -ConvertVector(Vector3D.CrossProduct(ConvertVector(unitVelocity), ConvertVector(unitLat)));
            Vector<double> unitDrag = -unitVelocity;

            Vector<double> forceAerodynamic = unitLift.Multiply(lift).Add(unitDrag.Multiply(drag));
            Vector<double> gravityForceN = new SparseVector(new[] {0, 0, m*g});
            
            Vector<double> force = forceAerodynamic.Add(gravityForceN); 
            
            Vector<double> mLat = unitLat.Multiply(adp*diameter*(CM + CMq*omegaLat));
            Vector<double> mSpin = new SparseVector(new [] {0, 0, CNr*omegaSpin}); // TODO: Check if missing element
            Vector<double> mvpInC = transformation.Multiply(mvp);
            Vector<double> mLatInC = transformation.Multiply(mLat);

            Vector<double> moment = mvpInC.Add(mLatInC).Add(mSpin);

            Vector<double> acceleration = force.Divide(m);


            double[] result = new double[12];

            result[0] = velocity[0];
            result[1] = velocity[1];
            result[2] = velocity[2];
            result[3] = acceleration[0];
            result[4] = acceleration[1];
            result[5] = acceleration[2];
            result[6] = st.PhiDot;
            result[7] = st.ThetaDot;
            result[8] = (moment[0] + Id*st.ThetaDot*st.PhiDot*st.SinTheta -
                         Ia*st.ThetaDot*(st.PhiDot*st.SinTheta + st.GammaDot) + Id*st.ThetaDot*st.PhiDot*st.SinTheta)/Id/
                        st.CosTheta;
            result[9] = (moment[1] + Ia*st.PhiDot*st.CosTheta*(st.PhiDot*st.SinTheta +st.GammaDot) - Id*st.PhiDot*st.PhiDot*st.CosTheta*st.SinTheta)/Id;
            result[10] = (moment[2] - Ia*(result[9]*st.SinTheta + st.ThetaDot*st.PhiDot*st.CosTheta))/Ia;
            result[11] = result[10];
            return result;
        }

        private static Vector3D ConvertVector(Vector<double> vector)
        {
            return new Vector3D(vector[0], vector[1], vector[2]);
        }

        private static Vector<double> ConvertVector(Vector3D vector)
        {
            return new SparseVector(new []{vector.X, vector.Y, vector.Z});
        }
    }
}
