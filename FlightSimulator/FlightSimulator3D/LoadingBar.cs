using System;
using System.Collections.Generic;
using System.Diagnostics;
using Mogre;

namespace Mogre.Demo.ExampleApplication
{
    public class LoadingBar
    {
        protected RenderWindow window;
	    protected Overlay loadOverlay;
	    protected float initProportion;
	    protected ushort numGroupsInit;
	    protected ushort numGroupsLoad;
	    protected float progressBarMaxSize;
	    protected float progressBarScriptSize;
	    protected float progressBarInc;
	    protected OverlayElement loadingBarElement;
	    protected OverlayElement loadingDescriptionElement;
	    protected OverlayElement loadingCommentElement;

        public void Start(RenderWindow window, ushort numGroupsInit, ushort numGroupsLoad, float initProportion) 
        {
            this.window = window;
            this.numGroupsInit = numGroupsInit;
            this.numGroupsLoad = numGroupsLoad;
            this.initProportion = initProportion;

            // We need to pre-initialise the 'Bootstrap' group so we can use
            // the basic contents in the loading screen
            ResourceGroupManager.Singleton.InitialiseResourceGroup("Bootstrap");

            OverlayManager omgr = OverlayManager.Singleton;
            loadOverlay = omgr.GetByName("Core/LoadOverlay");
            if (loadOverlay == null)
                System.Windows.Forms.MessageBox.Show("Cannot find loading overlay", "ExampleLoadingBar::start", System.Windows.Forms.MessageBoxButtons.OK, System.Windows.Forms.MessageBoxIcon.Error);
            loadOverlay.Show();

            // Save links to the bar and to the loading text, for updates as we go
            loadingBarElement = omgr.GetOverlayElement("Core/LoadPanel/Bar/Progress");
		    loadingCommentElement = omgr.GetOverlayElement("Core/LoadPanel/Comment");
		    loadingDescriptionElement = omgr.GetOverlayElement("Core/LoadPanel/Description");

            OverlayElement barContainer = omgr.GetOverlayElement("Core/LoadPanel/Bar");
		    progressBarMaxSize = barContainer.Width;
		    loadingBarElement.Width = 0F;
            
            ResourceGroupManager.Singleton.ResourceGroupScriptingStarted += ResourceGroupScriptingStarted;
            ResourceGroupManager.Singleton.ScriptParseStarted += ScriptParseStarted;
            ResourceGroupManager.Singleton.ScriptParseEnded += ScriptParseEnded;
            ResourceGroupManager.Singleton.ResourceGroupLoadStarted += ResourceGroupLoadStarted;
            ResourceGroupManager.Singleton.ResourceLoadStarted += ResourceLoadStarted;
            ResourceGroupManager.Singleton.WorldGeometryStageStarted += WorldGeometryStageStarted;
            ResourceGroupManager.Singleton.WorldGeometryStageEnded += WorldGeometryStageEnded;            
        }

        public void Finish()
        {
            // hide loading screen
            loadOverlay.Hide();

            ResourceGroupManager.Singleton.WorldGeometryStageEnded -= WorldGeometryStageEnded;
            ResourceGroupManager.Singleton.WorldGeometryStageStarted -= WorldGeometryStageStarted;           
            ResourceGroupManager.Singleton.ResourceLoadStarted -= ResourceLoadStarted;
            ResourceGroupManager.Singleton.ResourceGroupLoadStarted -= ResourceGroupLoadStarted;
            ResourceGroupManager.Singleton.ScriptParseEnded -= ScriptParseEnded;
            ResourceGroupManager.Singleton.ScriptParseStarted -= ScriptParseStarted;
            ResourceGroupManager.Singleton.ResourceGroupScriptingStarted -= ResourceGroupScriptingStarted;            
        }

        void ResourceGroupScriptingStarted(string groupName, uint scriptCount)
        {
            Debug.Assert(numGroupsInit > 0, "You stated you were not going to init any groups, but you did! Divide by zero would follow...");
            progressBarInc = progressBarMaxSize * initProportion / (float)scriptCount;
            progressBarInc /= numGroupsInit;
            loadingDescriptionElement.Caption = "Parsing scripts...";
            window.Update();
        }

        void ScriptParseStarted(string scriptName, out bool skipThisScript)
        {
            loadingCommentElement.Caption = scriptName;
            window.Update();

            skipThisScript = false;
        }

        void ScriptParseEnded(string scriptName, bool skipped)
        {
            loadingBarElement.Width += progressBarInc;
            window.Update();
        }

        void ResourceGroupLoadStarted(string groupName, uint resourceCount)
        {
            Debug.Assert(numGroupsLoad > 0, "You stated you were not going to load any groups, but you did! Divide by zero would follow...");
            progressBarInc = progressBarMaxSize * (1 - initProportion) / (float)resourceCount;
            progressBarInc /= numGroupsLoad;
            loadingDescriptionElement.Caption = "Loading resources...";
            window.Update();
        }

        void ResourceLoadStarted(ResourcePtr resource)
        {
            loadingCommentElement.Caption = resource.Name;
            window.Update();
        }

        void WorldGeometryStageStarted(string description)
        {
            loadingCommentElement.Caption = description;
            window.Update();
        }

        void WorldGeometryStageEnded()
        {
            loadingBarElement.Width += progressBarInc;
            window.Update();
        }
    }
}