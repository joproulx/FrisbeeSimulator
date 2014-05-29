using System;
using System.Collections.Generic;
using System.Text;
using System.Windows.Media.Media3D;
using FlightSimulator;

namespace Mogre.Demo.CameraTrack
{


    class Program
    {


        static void Main(string[] args)
        {
            try
            {
                CameraTrackApplication app = new CameraTrackApplication();
                app.Go();
            }
            catch (System.Runtime.InteropServices.SEHException)
            {
                // Check if it's an Ogre Exception
                if (OgreException.IsThrown)
                    ExampleApplication.Example.ShowOgreException();
                else
                    throw;
            }
        }
    }

    class CameraTrackApplication : Mogre.Demo.ExampleApplication.Example
    {
        //double[] m_y0Initials = new double[]
        //                          {
        //                              -9.03E-01,
        //                              -6.33E-01,
        //                              -9.13E-01,
        //                              1.34E+01,
        //                              -4.11E-01,
        //                              1.12E-03,
        //                              -7.11E-02,
        //                              2.11E-01,
        //                              -1.49E+01,
        //                              -1.48E+00,
        //                              5.43E+01,
        //                              5.03E+00
        //                          };

        //private double[] m_y0Initials = new double[]
        //                                    {
        //                                        -0.903,
        //                                        - 0.633,
        //                                        - 0.913,
        //                                        13.4	,
        //                                        -0.411,
        //                                            0.08864,
        //                                            1.00865,
        //                                            0.64820,
        //                                            -14.9,
        //                                            -1.48,
        //                                            54.3,
        //                                            5.7378651685393267
        //                                    };

        private double[] m_y0Initials = new double[]
                                            {
                                                -0.903,
                                                - 0.633,
                                                - 0.913,
                                                12.95056179775281,
                                                -0.23122471910112374	,
                                                    0.10212359550561803	,
                                                    -0.643033707865169	,
                                                    -0.025955056179775782,	
                                                    -14.9	,
                                                    -1.48	,
                                                    54.3,
                                                    4.6
                                            };


        private double[] m_y0;
        private double[,] m_y;
        AnimationState animState = null;

        private readonly Frisbee m_frisbee = new Frisbee();
        private SceneNode headNode;

        private float m_currentTime;
        private bool m_fStarted;

        public CameraTrackApplication()
        {
            m_y0 = new double[m_y0Initials.Length];

            double tfinal = 5; //% length of flight
            double nsteps = 292;// % number of time steps for data

            double span = tfinal / nsteps;

            double[] x = new double[(int)nsteps];
            for (int i = 0; i < nsteps; i++)
            {
                x[i] = span * i;
            }

            for (int i = 0; i < m_y0Initials.Length; i++)
            {
                m_y0[i] = m_y0Initials[i];
            }

            m_y = m_frisbee.Ode(m_y0, x); 


        }

        protected override void HandleInput(FrameEvent evt)
        {
            inputKeyboard.Capture();

            if (inputKeyboard.IsKeyDown(MOIS.KeyCode.KC_M))
            {
                m_fStarted = true;
            }

            if (inputKeyboard.IsKeyDown(MOIS.KeyCode.KC_N))
            {
                m_fStarted = false;
                m_currentTime = 0;
            }

            base.HandleInput(evt);
        }

        bool FrameStarted(FrameEvent evt)
        {
            //animState.AddTime(evt.timeSinceLastFrame);

            if (m_fStarted)
            {
                m_currentTime += evt.timeSinceLastFrame;
                int index = GetIndex(evt);

                headNode.SetPosition(Convert(m_y[index, 1]), Convert(m_y[index, 2]), Convert(m_y[index, 3]));

                float f = 0f;
                if (m_currentTime <= 10)
                {
                    f = -m_currentTime/10F*Math.PI;
                }
                else
                {
                    int test = 0;
                }
                
                //var quaternion = new Quaternion(GetRotation((float)m_y[index, 7], (float)m_y[index, 8] + Math.PI / 2F, 0));
                var quaternion = new Quaternion(GetRotation(0, 0, f));
                headNode.ResetOrientation();
                headNode.Rotate(new Vector3(0, 0, 1F), Math.HALF_PI);
                headNode.Rotate(new Vector3(0, -1F, 0), (float)m_y[index, 8]);
                headNode.Rotate(new Vector3(0, 0, 1F), (float)m_y[index, 7]);
                
                //headNode.SetOrientation(quaternion.w, quaternion.x, quaternion.y, quaternion.z);
            
            }
            return true;

        }


        Matrix3 GetRotation(float roll, float yaw, float pitch)
        {
            // g = gamma = yaw
            // t = theta = pitch
            // p = phi   = roll
            //       --                                      --              
            //      |                                                                        cos(g)*sin(pitch)*cos(roll)+sin(yaw)*sin(roll) | 
            //      |                                                                                       |
            //  T = | sin(yaw)*cos(pitch), sin(yaw)*sin(pitch)*sin(roll)-cos(yaw)*cos(roll), sin(yaw)*sin(pitch)*cos(roll)+cos(yaw)*sin(roll) |
            //      |                                                                                       |
            //      |    -sin(pitch)   ,          cos(pitch)*sin(roll)            ,          cos(pitch)*cos(roll)             |                
            //       --                                                                                   --
            

            return new Matrix3(new [,]
                                   {
                                       { Math.Cos(yaw)*Math.Cos(pitch), Math.Cos(yaw)*Math.Sin(pitch) * Math.Sin(roll) - Math.Sin(yaw) * Math.Cos(roll), Math.Cos(yaw)* Math.Sin(pitch) * Math.Cos(roll) + Math.Sin(yaw)*Math.Sin(roll) },
                                       {Math.Sin(yaw) * Math.Cos(pitch), Math.Sin(yaw)*Math.Sin(pitch)*Math.Sin(roll)-Math.Cos(yaw)*Math.Cos(roll), Math.Sin(yaw)*Math.Sin(pitch)*Math.Cos(roll) +Math.Cos(yaw)* Math.Sin(roll)},
                                       { -Math.Sin(pitch), Math.Cos(pitch) * Math.Sin(roll), Math.Cos(pitch) * Math.Cos(roll) }
                                   });

        }

        private float Convert(double p)
        {
            return (float)p*10F;
        }
        
        private int GetIndex(FrameEvent evt)
        {
            int i = 0;
            for (int j = 0; j < m_y.GetLength(0); j++)
            {
                if (m_currentTime < m_y[j, 0] || i == (m_y.GetLength(0) - 1))
                {
                    break;
                }
                i++;
            }

            return i;
        }

        public override void CreateFrameListener()
        {
            base.CreateFrameListener();
            Root.Singleton.FrameStarted += FrameStarted;
        }

        // Scene creation
        public override void CreateScene()
        {
            // Set ambient light
            sceneMgr.AmbientLight = new ColourValue(0.2F, 0.2F, 0.2F);

            // Create a skydome
            sceneMgr.SetSkyDome(true, "Examples/CloudySky", 5, 8);

            // Create a light
            Light l = sceneMgr.CreateLight("MainLight");

            // Accept default settings: point light, white diffuse, just set position
            // NB I could attach the light to a SceneNode if I wanted it to move automatically with
            //  other objects, but I don't
            l.Position = new Vector3(20F, 80F, 50F);

            // Define a floor plane mesh
            Plane p;
            p.normal = Vector3.UNIT_Y;
            p.d = 200;
            MeshManager.Singleton.CreatePlane("FloorPlane", ResourceGroupManager.DEFAULT_RESOURCE_GROUP_NAME, p, 200000F, 20000F, 20, 20, true, 1, 50F, 50F, Vector3.UNIT_Z);

            Entity ent;
            // Create an entity (the floor)
            ent = sceneMgr.CreateEntity("floor", "FloorPlane");
            ent.SetMaterialName("Examples/RustySteel");

            // Attach to child of root node, better for culling (otherwise bounds are the combination of the 2)
            sceneMgr.RootSceneNode.CreateChildSceneNode().AttachObject(ent);

            // Add a head, give it it's own node
            headNode = sceneMgr.RootSceneNode.CreateChildSceneNode();
            ent = sceneMgr.CreateEntity("head", "cylinder.mesh");
            //ent.SetMaterialName("materials/texturesRustySteel");

            //       y
            //       ^
            //       |       
            // x <---x z


           
            
            headNode.AttachObject(ent);
            headNode.Scale(2,10,10);
            //headNode.Rotate(new Vector3(0,0,1), new Radian(Math.PI / 2F));
            headNode.Rotate(new Vector3(0, 0, 1F), Math.HALF_PI);
                
            
            // Make sure the camera track this node
            //camera.SetAutoTracking(true, headNode);

            // Create the camera node & attach camera
            SceneNode camNode = sceneMgr.RootSceneNode.CreateChildSceneNode();
            camNode.AttachObject(camera);

            // set up spline animation of node
            //Animation anim = sceneMgr.CreateAnimation("CameraTrack", 10F);

            // Spline it for nice curves

            //anim.SetInterpolationMode(Animation.InterpolationMode.IM_SPLINE);
            //// Create a track to animate the camera's node
            //NodeAnimationTrack track = anim.CreateNodeTrack(0, camNode);

            //// Setup keyframes
            //TransformKeyFrame key = track.CreateNodeKeyFrame(0F); // startposition
            //key = track.CreateNodeKeyFrame(2.5F);
            //key.Translate = new Vector3(500F, 500F, -1000F);
            //key = track.CreateNodeKeyFrame(5F);
            //key.Translate = new Vector3(-1500F, 1000F, -600F);
            //key = track.CreateNodeKeyFrame(7.5F);
            //key.Translate = new Vector3(0F, 100F, 0F);
            //key = track.CreateNodeKeyFrame(10F);
            //key.Translate = new Vector3(0F, 0F, 0F);

            //// Create a new animation state to track this
            //animState = sceneMgr.CreateAnimationState("CameraTrack");
            //animState.Enabled = true;

            // Put in a bit of fog for the hell of it        
            //sceneMgr.SetFog(FogMode.FOG_EXP, ColourValue.White, 0.0002F);
        }
    }
}