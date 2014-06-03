using System;
using System.Collections.Generic;
using System.Text;

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
        AnimationState animState = null;

        bool FrameStarted(FrameEvent evt)
        {
            animState.AddTime(evt.timeSinceLastFrame);
            return true;
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
            SceneNode headNode = sceneMgr.RootSceneNode.CreateChildSceneNode();
            ent = sceneMgr.CreateEntity("head", "ogrehead.mesh");
            
            headNode.AttachObject(ent);

            // Make sure the camera track this node
            camera.SetAutoTracking(true, headNode);

            // Create the camera node & attach camera
            SceneNode camNode = sceneMgr.RootSceneNode.CreateChildSceneNode();
            camNode.AttachObject(camera);

            // set up spline animation of node
            Animation anim = sceneMgr.CreateAnimation("CameraTrack", 10F);
    
            // Spline it for nice curves
    
            anim.SetInterpolationMode(Animation.InterpolationMode.IM_SPLINE);
            // Create a track to animate the camera's node
            NodeAnimationTrack track = anim.CreateNodeTrack(0, camNode);
    
            // Setup keyframes
            TransformKeyFrame key = track.CreateNodeKeyFrame(0F); // startposition
            key = track.CreateNodeKeyFrame(2.5F);
            key.Translate = new Vector3(500F, 500F, -1000F);
            key = track.CreateNodeKeyFrame(5F);
            key.Translate = new Vector3(-1500F, 1000F, -600F);
            key = track.CreateNodeKeyFrame(7.5F);
            key.Translate = new Vector3(0F, 100F, 0F);
            key = track.CreateNodeKeyFrame(10F);
            key.Translate = new Vector3(0F, 0F, 0F);

            // Create a new animation state to track this
            animState = sceneMgr.CreateAnimationState("CameraTrack");
            animState.Enabled = true;        
    
            // Put in a bit of fog for the hell of it        
            sceneMgr.SetFog(FogMode.FOG_EXP, ColourValue.White, 0.0002F);            
        }
    }
}