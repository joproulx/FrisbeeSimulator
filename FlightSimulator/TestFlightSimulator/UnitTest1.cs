using System;
using System.Text;
using System.Collections.Generic;
using System.Linq;
using FlightSimulator;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using MathNet.Numerics.Interpolation.Algorithms;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics.LinearAlgebra.Double.Factorization;
using MathNet.Numerics.LinearAlgebra.Generic;

namespace TestFlightSimulator
{
    [TestClass]
    public class UnitTest1
    {
        [TestMethod]
        public void TestVelocityOnPlane()
        {
            Frisbee.SimulationState st = new Frisbee.SimulationState
                                             {
                                                 VX = 1,
                                                 VY = 1,
                                                 Theta = Math.PI / 4
                                             };


            Matrix<double> transformation =
                new SparseMatrix(new [,]
                                     {
                                         {st.CosTheta, st.SinTheta*st.SinPhi, -st.SinTheta*st.CosPhi}, 
                                         {0, st.CosPhi, st.SinPhi}, 
                                         {st.SinTheta, -st.CosTheta*st.SinPhi, st.CosTheta*st.CosPhi}
                                     });


            SparseVector c3 = new SparseVector(transformation.Row(2));

            SparseVector velocity = new SparseVector(new[] { st.VX, st.VY, st.VZ });
            double velocityMagnitude = velocity.Norm(2);

            double velocityC3 = velocity.DotProduct(c3);

            Vector<double> vp = velocity.Subtract(c3.Multiply(velocityC3));
            double vpMagnitude = vp.Norm(2);


        }
    }
}
