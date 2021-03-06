import sys,os,io

cur = os.path.dirname( os.path.realpath( __file__ ) )

sys.path.append(cur+'/../UtilityLib')

from moduleXML import loadNetworkFromXML 
from moduleStartUp import parseOptions
from modulePickle import loadSolutionDataFile
from modulePickle import loadExternalDataSet

from processing import linearWaveSplitting
from processing import nonLinearWaveSplitting
from processing import minMaxFunction

sys.path.append(cur+'/../NetworkLib')
from classVascularNetwork import VascularNetwork 
from classVessel import Vessel 

from class3dLUT import WindowLUT
from class3dLUT import LUT
from class3dControlGUI import ControlWindow

import numpy as np
from math  import *

import copy 
import time
import gc

from OpenGL.GLUT import *
from OpenGL.GL import *
from OpenGL.GLU import * 

import pyglet
from pyglet.gl import *
from pyglet import *
from pyglet.window import mouse, key
from pyglet.window.key import MOTION_UP,MOTION_DOWN

import matplotlib.animation as animation
from matplotlib._png import read_png
from pylab import *    

__author__ = "Vinzenz Gregor Eck"
__version__ = "0.3"
  
class Visualisation3DGUI(pyglet.window.Window):
    '''
    This class defines the GuI inclusive call buttons
    (saparete GUI class to Gui-user Class
    which  is then wrapped into Visualisation3D )
    '''
    def __init__(self):

        ## openGL context configuration
        stencil_size    = 8
        depth_size      = 16
        samples         = 4
        sample_buffers  = 1
        
        platform = pyglet.window.get_platform()
        display = platform.get_default_display()
        screen = display.get_default_screen()
        
        templateHigh = pyglet.gl.Config(sample_buffers=sample_buffers,
                                    samples = samples,
                                    double_buffer=True,
                                    depth_size = depth_size,
                                    stencil_size = stencil_size)
        
        templateLow = pyglet.gl.Config(sample_buffers=0,
                                    double_buffer=False,
                                    depth_size = 16)
        
        try:
            config = screen.get_best_config(templateHigh)
        except pyglet.window.NoSuchConfigException:
            try:
                config = screen.get_best_config(templateLow)
                print "OpenGlconfig", config
            except pyglet.window.NoSuchConfigException:
                print " OpenGL: simple context configuration is applied due to installed hardware"
                template = gl.Config()
                config = screen.get_best_config(template)
                print "OpenGlconfig", config
        
        context = config.create_context(None)
        
        super(Visualisation3DGUI, self).__init__(resizable = True, context = context)
                        
        self.set_maximum_size(screen.width, screen.height)
        self.set_minimum_size(self.width, self.height)
                     
        self.keys = key.KeyStateHandler()
        self.push_handlers(self.keys)
        
        
        self.backgroundColorIterator = 0
        self.backgroundColors = [[1.0,1.0,1.0,1],
                                 [0.95,0.95,0.95,1],
                                 [0.85,0.85,0.85,1],
                                 [0.75,0.75,0.75,1]]
        self.backgroundColor = [0.75,0.75,0.75,1]
                        
        pyglet.gl.glClearColor(*self.backgroundColor)
        
        
        self.batch = pyglet.graphics.Batch()
        self.beta  = (pi*5.)/180.
        self.alpha = (pi*5.)/180.
        
        self.dalpha = (pi*5.)/180.
        
        
        self.cameraEyesInit   = np.array([4.,0.,0.])
        self.cameraLookAtInit = np.array([0.,0.,0.])
        self.cameraEyes = np.copy(self.cameraEyesInit)
        self.cameraLookAt = np.copy(self.cameraLookAtInit)
        self.normal = np.array([0.,0.,0.])
        #self.viewUp  = np.array([-0.70710678 , 0., 0.70710678]) # 45 degree
        self.viewUp  = np.array([0 , 0., 1])
        self.axisX = np.array([1.,0.,0])
        self.axisY = np.array([0.,1.,0])
        self.axisZ = np.array([0.,0.,1])
        
        self.viewZoomingSensitivity   = 0.5
        self.viewTranslateSensitivityX = 0.001
        self.viewTranslateSensitivityY = 0.0005
        
        self.viewRotationSensitivityZ = self.dalpha ###
        
        # simulation runtime variables
        self.simulationRunning = False
        self.clock = pyglet.clock.get_default()
        self.updateTime = 1./60.
                
        # One-time GL setup
        glColor3f(1, 0, 0)
        glEnable(GL_DEPTH_TEST)
        ## uncomment this line for viewing only front faces
        #glEnable(GL_CULL_FACE)
        
        # Uncomment this line for a wireframe view
        #glPolygonMode(GL_FRONT_AND_BACK, GL_LINE)
        
        # Simple light setup.  On Windows GL_LIGHT0 is enabled by default,
        # but this is not the case on Linux or Mac, so remember to always 
        # include it.
        glEnable(GL_LIGHTING)
        glEnable(GL_LIGHT0)
        glEnable(GL_LIGHT1)
        
        # Define a simple function to create ctypes arrays of floats:
        def vec(*args):
            return (GLfloat * len(args))(*args)
        
        glLightfv(GL_LIGHT0, GL_POSITION, vec(.5, .5, 1, 0))
        glLightfv(GL_LIGHT0, GL_SPECULAR, vec(.5, .5, 1, 1))
        glLightfv(GL_LIGHT0, GL_DIFFUSE, vec(1, 1, 1, 1))
        glLightfv(GL_LIGHT1, GL_POSITION, vec(1, 0, .5, 0))
        glLightfv(GL_LIGHT1, GL_DIFFUSE, vec(.5, .5, .5, 1))
        glLightfv(GL_LIGHT1, GL_SPECULAR, vec(1, 1, 1, 1))
        
        
        ###
        # control GUI window
        self.controlWindow = ControlWindow(self)
        
        # register button events
        # simulation speed
        self.controlWindow.onPressRealTimeSim       = self.visualisationSpeedRealTime
        self.controlWindow.onPressHalfTimeSim       = self.visualisationSpeedHalfTime
        self.controlWindow.onChangeSpeedSlider      = self.setVisualisationSpeed
        # view
        self.controlWindow.onPressViewXY            = self.setViewXY
        self.controlWindow.onPressViewYZ            = self.setViewYZ
        self.controlWindow.onPressViewXZ            = self.setViewXZ
        self.controlWindow.onPressViewXYZ           = self.setViewXYZ
        # play pause / reset
        self.controlWindow.onPressPlayPause         = self.startVisualisation
        self.controlWindow.onPressReset             = self.resetVisualisationAndPause
        # create movie
        self.controlWindow.onPressMovie             = self.recordMovie
        self.controlWindow.onPressViewPhoto         = self.saveScreenShot
        # LUT    
        self.controlWindow.onPressLookUpTable       = self.changeColorTable    
        self.controlWindow.onPressChangeQuantity    = self.changeQuantityLUT
        # wave split
        self.controlWindow.onPressSplit             = self.enableWaveSplitLookUp
        self.controlWindow.onPressSplitExp          = self.enableLeftRightColoring
        # vessel area/wall movement
        self.controlWindow.onPressWallMovement      = self.enableWallMovement
        self.controlWindow.onChangeWallSlider       = self.setAreaFactor
        self.controlWindow.onPressBackgroundColor   = self.changeBackgroundColor
        
    def on_resize(self,width, height):
        self.switch_to()
        glViewport(0, 0, width, height)
        self.camera()
        return pyglet.event.EVENT_HANDLED
                
    def on_draw(self):
        self.switch_to()
        self.clear() # == glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT )
        pyglet.gl.glClearColor(*self.backgroundColor)
        self.camera()
        #self.drawPlane()
        #self.drawCoordinateSystem()
        glEnable(GL_COLOR_MATERIAL)
        self.batch.draw()
        return pyglet.event.EVENT_HANDLED
    
    def resizeWindow(self,dx,dy):
        print "new size", self.width+dx, self.height-dy
        self.set_size(self.width+dx, self.height-dy) 
        
    def drawPlane(self):
        self.planeZ = -0.25
        xWidth = 1
        yWidth = 1
        glBegin(GL_QUADS)
        glColor3f(0.8,0.8,0.8)
        glVertex3f(xWidth,yWidth,self.planeZ)
        glNormal3f(0, 0, abs(self.planeZ/self.planeZ))
        glVertex3f(xWidth, -yWidth,self.planeZ)
        glNormal3f(0, 0, abs(self.planeZ/self.planeZ))
        glVertex3f(-xWidth, -yWidth,self.planeZ)
        glNormal3f(0, 0, abs(self.planeZ/self.planeZ))
        glVertex3f(-xWidth, yWidth,self.planeZ)
        glNormal3f(0, 0, abs(self.planeZ/self.planeZ))
        glEnd()
                
    def drawCoordinateSystem(self):
        glClear(GL_COLOR_BUFFER_BIT)
        glBegin( GL_LINES )
        glColor3f( 1, 0, 0 )
        glVertex3f( 0, 0, 0 )   
        glVertex3f( 0.1, 0, 0 )
        glEnd( )
        
        glBegin( GL_LINES )
        glColor3f(0,1,0)
        glVertex3f( 0, 0, 0 )   
        glVertex3f( 0, 0.1, 0 )
        glEnd( )
     
        glBegin( GL_LINES )
        glColor3f(0,0,1)
        glVertex3f( 0, 0, 0 )   
        glVertex3f( 0, 0, 0.1 )
        glEnd( )
        glFlush()
    
    def camera(self):
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity() 
        gluPerspective(5,1,1,1000)
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()
        gluLookAt(self.cameraEyes[0],self.cameraEyes[1],self.cameraEyes[2],
                  self.cameraLookAt[0],self.cameraLookAt[1],self.cameraLookAt[2],
                  self.viewUp[0],self.viewUp[1],self.viewUp[2])
       
    def cameraTranslateZ(self,dx,dy):
        """
        click left mouse and move 3D vessel around
        """
        kz = -dy*self.viewTranslateSensitivityY
                
        self.cameraEyes   += self.axisZ*kz
        self.cameraLookAt += self.axisZ*kz
             
    def cameraTranslateXY(self,dx,dy):
        """
        click right mouse and move vessel up and down
        """
        ky = dy*self.viewTranslateSensitivityY
        kx = dx*self.viewTranslateSensitivityX

        if self.cameraEyes[2]<0: ky = -ky
        if (self.cameraEyes[0:2] == self.cameraLookAt[0:2]).all():
            translation = - ky* self.axisX + kx*self.axisY

        else:
            eyesLookAtVector = (self.cameraEyes-self.cameraLookAt)*(self.axisX+self.axisY)
            normal =  np.cross(eyesLookAtVector,self.axisZ)
            normal = normal/np.linalg.norm(normal)
            translation = kx*normal + ky*eyesLookAtVector
              
        self.cameraEyes   = self.cameraEyes   + translation
        self.cameraLookAt = self.cameraLookAt + translation

        
    def zoomInAndOut(self,dx,dy):
        """
        click middle mouse and zoom in and out
        """        
        eyesLookAtVector = self.cameraLookAt-self.cameraEyes
        normLE = np.linalg.norm(eyesLookAtVector)
        
        ## zooming out
        if dx < dy :
            self.cameraEyes -= self.viewZoomingSensitivity*eyesLookAtVector/normLE
        ## zooming in
        elif dx > dy:
            if normLE > 2:
                self.cameraEyes += self.viewZoomingSensitivity*eyesLookAtVector/normLE

    
    def rotate(self,dx,dy):
        """
        click right and left mouse to rotate around
        the z axis and the axis normal to the view plane
        """  
        eyesLookAtVector = self.cameraEyes-self.cameraLookAt # EL
        #print 'eyesLookAtVector', eyesLookAtVector
        eyesLookAtVectorProjection = eyesLookAtVector*np.array([1.,1.,0.]) #ELxy
          
        n = np.dot(eyesLookAtVectorProjection,eyesLookAtVector)
        normEyesLookAtVector = np.linalg.norm(eyesLookAtVector)
        normEyesLookAtVectorProjection = np.linalg.norm(eyesLookAtVectorProjection)
          
        normProduct = normEyesLookAtVector*normEyesLookAtVectorProjection 
        if normProduct == 0: normProduct = 1
        
        cosAlpha = (n/normProduct)
        if cosAlpha > 1 : cosAlpha = 1.0
        elif cosAlpha < -1 : cosAlpha = -1.0 
                
        alpha =  np.int(np.floor(np.arccos(cosAlpha)/np.pi*180.))
        
        degreeFactor = 180/pi
        radianFactor = pi/180
          
                  
        if alpha < 1.:   # <---------------------------------------------------------CASE 1
            self.viewUp= np.array([0,0,1])
            if abs(dx) > 2:  self.cameraRotateAroundZ(dx)
            if abs(dy) > 1: self.cameraRotateGear(dy)
            
        elif alpha >= 89.:#<--------------------------------------------------------CASE 2
            if abs(dy) > 1: self.cameraRotateGear(dy)
            if abs(dx) > 2: self.cameraRotateAroundZ(dx, rotateViewUp = True)
            
        else: #<----------------------------------------------------------------------CASE 3
            # rotate
            n = np.cross(eyesLookAtVectorProjection,eyesLookAtVector)
            viewUp = np.cross(n,eyesLookAtVector)
            viewUpNorm = np.linalg.norm(viewUp)
            if self.cameraEyes[2] > self.cameraLookAt[2]: self.viewUp  = viewUp/viewUpNorm
            else: self.viewUp  = - viewUp/viewUpNorm
            
            if abs(dx) > 2: self.cameraRotateAroundZ(dx)
            if abs(dy) > 1: self.cameraRotateGear(dy)
                 
    def cameraRotateAroundZ(self, dx, rotateViewUp = False):
        signDx = -dx/abs(dx)
        Rz = np.array( [[cos(self.alpha),        signDx*sin(self.alpha), 0.],
                        [-signDx*sin(self.alpha), cos(self.alpha),        0.],
                        [0.,                     0.,                     1.]])
        if rotateViewUp == False:
            self.cameraEyes = np.dot(self.cameraEyes,Rz)
        else:        
            self.viewUp = np.dot(self.viewUp,Rz)
           
    def cameraRotateGear(self,dy):     
                ## if mouse moves horizontaly
                        
                s = dy/abs(dy)
    
                #E = np.array( [self.EX,self.EY,0] )
                E = self.cameraEyes*(self.axisX+self.axisY)
                    
                self.normal  = np.cross(self.axisZ,E) # normal to plane (z,E)
                lengthN = np.linalg.norm(self.normal)
                
                ux, uy, uz = -self.normal[0]/lengthN, -self.normal[1]/lengthN, -self.normal[2]/lengthN # normalize normal coordinates
                
                #if n == 0:
                    
                
                A11 = cos(self.beta)+ux*ux*(1-cos(self.beta))
                A12 = ux*uy*(1-cos(self.beta))-uz*sin(self.beta)
                A13 = ux*uz*(1-cos(self.beta))+uy*sin(self.beta)
                    
                A21 = uy*ux*(1-cos(self.beta))+uz*sin(self.beta)
                A22 = cos(self.beta)+uy*uy*(1-cos(self.beta))
                A23 = uy*uz*(1-cos(self.beta))-ux*sin(self.beta)
                    
                A31 = uz*ux*(1-cos(self.beta))-uy*sin(self.beta)
                A32 = uz*uy*(1-cos(self.beta))+ux*sin(self.beta)
                A33 = cos(self.beta)+uz*uz*(1-cos(self.beta))
                    
                Rn = np.array( [[A11, A12, A13],
                                [A21, A22, A23], 
                                [A31, A32, A33]] )
                              
                if s < 0:   self.cameraEyes = np.dot(self.cameraEyes,Rn) # upwards rotation
                elif s > 0: self.cameraEyes = np.dot(Rn,self.cameraEyes.T).T #downwards
                        
# 
            
    def on_mouse_drag(self,x,y,dx,dy,buttons,modifiers):
        
        if buttons == 1: # left mouse button
            if self.keys[key.LCTRL] == False:
                self.cameraTranslateXY(dx,dy)
            else:
                self.resizeWindow(dx, dy)
            
        elif buttons == 4: # right mouse button
            self.cameraTranslateZ(dx,-dy)
            
        elif buttons == 5: # left plus right
            self.rotate(dx,dy)
                   
        elif buttons == 2: # middle button to zoo
            self.zoomInAndOut(dx,dy)
            
    def on_mouse_scroll(self, x, y, scroll_x, scroll_y):
        self.zoomInAndOut(0,-scroll_y)
        
    def on_key_press(self,symbol, modifiers):
        
        
        #### start stop simulation visualisation
        if symbol == key.P:
            self.startVisualisation()
        elif symbol == key.O:
            self.resetVisualisationAndPause()
        
        
                            
        #### visualisation speed
        elif symbol == key.H:
            self.visualisationSpeedHalfTime()
        elif symbol == key.R:
            self.visualisationSpeedRealTime()
        elif symbol == key.PLUS:
            self.visualisationSpeedUp()
        elif symbol == key.MINUS:
            self.visualisationSpeedDown()
            
        #### look up table stuff
        elif symbol == key.L:
            self.changeColorTable()
        
        elif symbol == key.W:
            self.enableWaveSplitLookUp()
            
        elif symbol == key.E:
            self.enableLeftRightColoring()
            
        elif symbol == key.Q:
            self.changeQuantityLUT()
            
        elif symbol == key.A:
            self.enableWallMovement()
            
        ### create movie:
        elif symbol == key.M:
            self.recordMovie()
            
        elif symbol == key.S:
            self.saveScreenShot()
        
        elif symbol == pyglet.window.key.ESCAPE:
            self.on_close()
                    
        elif symbol == key.B:
            self.changeBackgroundColor()
                   
        ## view 
        elif symbol == key._1 :
            self.setViewXY()
        elif symbol == key._2 :
            self.setViewYZ()
        elif symbol == key._3 :
            self.setViewXZ()
        elif symbol == key._4 :
            self.setViewXYZ()()     
          
        elif symbol == key.NUM_8:
            # rotate up wards
            self.rotate(0,2)
        elif symbol == key.NUM_2:
            # rotate downwards
            self.rotate(0,-2)
        elif symbol == key.NUM_6:
            # rotate right
            self.rotate(4,0)
        elif symbol == key.NUM_4:
            # rotate left
            self.rotate(-4,0)
            
            
    def on_close(self):
        self.cLUT.windowLUT.switch_to()
        pyglet.app.exit()
        self.controlWindow.switch_to()
        pyglet.app.exit()
        self.switch_to()
        pyglet.app.exit()
      
    ### callback functions
        
    def startVisualisation(self):
        if self.simulationRunning: #stop it
                self.clock.unschedule(self.updateVisualisation)
                self.simulationRunning = False
        else: #start it
            self.clock.schedule_interval(self.updateVisualisation,self.updateTime) #
            self.simulationRunning = True
    
    def resetVisualisationAndPause(self):
        if self.simulationRunning : self.startVisualisation()
        self.timeStepCurrent = 0
        self.updateVisualisation(self.updateTime)
        self.timeStepCurrent = 0
        self.switch_to()
        self.on_draw()
        
        
    def updateVisualisation(self,dt):
        pass   
    
    ## visualisation speed
    def visualisationSpeedHalfTime(self):
        pass
    
    def visualisationSpeedRealTime(self):
        pass
        
    def visualisationSpeedUp(self):
        pass
    
    def visualisationSpeedDown(self):
        pass
    
    def setVisualisationSpeed(self, value):
        pass
    
    def changeColorTable(self):       
        pass
    
    ### area / vessel wall
    def enableWallMovement(self):
        pass
    
    def setAreaFactor(self):
        pass
    
    ### wave split
    def enableWaveSplitLookUp(self):
        pass
     
    def enableLeftRightColoring(self):
        pass
     
    def changeQuantityLUT(self):
        pass
    
    ### record save movie screenshot
    def recordMovie(self):
        pass

    def saveScreenShot(self):
        pass
        
    ###
    def changeBackgroundColor(self):
        tableIterator = range(len(self.backgroundColors))
        tableIterator.pop(0)
        tableIterator.append(0)      
        self.backgroundColor = self.backgroundColors[self.backgroundColorIterator]
        self.backgroundColorIterator = tableIterator[self.backgroundColorIterator]
     
    ### view
    def setViewXY(self):
        self.cameraEyes   = np.array([0.,0.,4.])
        self.cameraLookAt = np.array([0.,0.,0.])
        self.viewUp = np.array([1.,0.,0.])
            
    def setViewYZ(self):
        self.cameraEyes   = np.array([4.,0.,0.])
        self.cameraLookAt = np.array([0.,0.,0.])
        self.viewUp = np.array([0.,0.,1.])

    def setViewXZ(self):
        self.cameraEyes   = np.array([0.,-4.,0.])
        self.cameraLookAt = np.array([0.,0.,0.])
        self.viewUp = np.array([0.,0.,1.])
    
    def setViewXYZ(self):
        self.cameraEyes = np.copy(self.cameraEyesInit)
        self.cameraLookAt = np.copy(self.cameraLookAtInit)
        self.viewUp = np.array([0.,0.,1.])
        
        
class Visualisation3D(Visualisation3DGUI):
    '''
    This class defines the 3D visualisation
    it is an inheritance of Visualisation3DGui which defines the GUI
    '''
    def __init__(self, networkName = None, dataSetNumber = None, connect = False):
        
        Visualisation3DGUI.__init__(self)
        
        # load functions
        self.networkName = networkName
        self.dataSetNumber = dataSetNumber
        self.connect = connect
               
        # vascularNewtor Instance
        self.vascularNetwork = None
        
        # dictionary with 3D vessels
        self.vessels3D = {}
        self.nPointsPerCircle = 8  # points for the circle at one gridpoint should bee around 8-16-32
        
        # vessel area/wall 
        self.areaFactor         = [2]
        self.areaFactorMaximum  = 200
        self.wallMovement       = [True]
            
        # look up table 
        self.pressureMaxima      = None
        self.pressureMaximaSplit = None
        self.cLUT = None # color look up table class
        self.quantityLUT = ['Pressure'] 
                
        # enable waveSplit LUT
        self.waveSplit = [False] # us an array as it is behaving like a pointer
        
        # simulation timing and update 
        self.updateTime = 1./60.
        self.totalTime = 0
        self.timeStepsTotal     = 50 # get this from solution length e.g. len(area)
        self.timeStepCurrent    = 0
        self.timeStepIncrement  = 1 # time step increment from [1, ... self.totalTimeSteps-1 ]
                
        ## create movie
        self.recordMovieData = False
        self.movieCount = 0
        self.movieCount = 0
        
        self.templateMovDataDirectory = ''.join([cur,'/.movieTemplateData'])
        self.movieSaveDirectory = ''.join([cur+'/Movies'])
        try: os.mkdir(self.movieSaveDirectory)
        except: pass
        self.movieNumber = 0
        self.movieFileType = '.mp4'
        
        ## create screen shot
        self.screenShotNumber = 0 
        self.screenShotDataDirectory = ''.join([cur,'/screenShots'])
                
        self.initialzeVisualisation()
        
    
    #--- functions for initialisation
    def initialzeVisualisation(self):
        '''
        
        '''        
        # to be changed
        self.loadVascularNetwork()
        # find min max pressures and apply it to LUT
        self.findMinMaxPressure()
        # create LUT      
        self.cLUT = LUT(self.quantitiyMaxima) 
        # creat 3d vessel representations
        
        self.createVessel3D()
        start = time.clock()
        # apply number total timesteps
        self.timeStepsTotal = len(self.vascularNetwork.vessels[self.vascularNetwork.root].Psol)
        
        # evaluate maximum area factor
        self.findMaximumAreaFactor()
        
        # evaluate visualisation time step
        self.evaluateUpdateRate()
        
        # small stuff
        self.set_caption(self.vascularNetwork.name)
        
        self.networkName = self.vascularNetwork.name   
        # evaluate movie count from already saved movies:
        for dirName, dirNames, fileNames in os.walk(self.movieSaveDirectory):
            for fileName in fileNames:
                if self.networkName in fileName:
                    try:
                        number = int(fileName.split('_')[-1].split(self.movieFileType)[0])
                        self.movieNumber = max(self.movieNumber,1 + number)
                    except: pass
        # evaluate the screenshot count from already saved screenshots
        for dirName, dirNames, fileNames in os.walk(self.screenShotDataDirectory):
            for fileName in fileNames:
                if self.networkName in fileName:
                    try: 
                        number = int(fileName.split('_')[-1].split('.png')[0])
                        self.screenShotNumber = max(self.screenShotNumber,1 + number)
                    except: pass
        
    def loadVascularNetwork(self):

        '''
        load vascularNetwork inclusiv all solution data
        from solution data set
        
        uses local variables networkName and dataSetNumber
        saves vascularNetwork instance in self.vascularNetwork
        '''
        if self.networkName and self.dataSetNumber:
            vascularNetwork, solutionDataSets, simulationCaseDescriptions = loadSolutionDataFile(self.networkName, self.dataSetNumber) 
            self.vascularNetwork = vascularNetwork    
        
    def createVessel3D(self):

        '''
        method to create a Vessel3D instance for each vessel in the current
        vascularNetwork and stores it in self.vessels3D = {vesselId: vessel3D}
        '''
        if self.vascularNetwork is not None:  
            for vesselId,vessel in self.vascularNetwork.vessels.iteritems():
                self.vessels3D[vesselId] = Vessel3D( vessel,
                                                     self.batch,
                                                     self.nPointsPerCircle,
                                                     self.wallMovement,
                                                     self.areaFactor,
                                                     self.cLUT,
                                                     self.quantityLUT,
                                                     self.waveSplit)
       
    #--- functions for update                 
    def updateVisualisation(self, dt):

        ''' 
        method to update the visualisation of the 
        network depending on the timeStep
        this method is call continously if the "play button is clicked"
        '''
        timeStepCurrent = self.timeStepCurrent     
        ## update all visualisations of all vessel3D instances   
        # update vessel Positions, 3dgeometry and colour of vertices
        for vesselId in self.vascularNetwork.treeTraverseList:
            self.vessels3D[vesselId].update3Dposition(timeStepCurrent,None, None)
        
        # update the time step for next call
        self.timeStepCurrent = self.timeStepCurrent+self.timeStepIncrement
        # restart simulation if end reached
        if self.timeStepCurrent >= self.timeStepsTotal:
            self.timeStepCurrent = 0
            self.timeElapsed = 0
        
        # record movie data
        if self.recordMovieData == True:
            ## save screenshots
            self.switch_to()
            pyglet.image.get_buffer_manager().get_color_buffer().save(''.join([self.templateMovDataDirectory,'/.screenShot',str(self.movieCount),'.png']))
            self.movieCount = self.movieCount + 1 
        
    def recordMovie(self):
        '''
        start end end recording movies
        '''
        # start recording
        if self.recordMovieData == False:
                       
            try: os.mkdir(self.templateMovDataDirectory)
            except: pass
            
            self.movieCount = 0
            self.recordMovieData = True
            
        else: # create movie out of recorded data 
            self.recordMovieData = False
            self.createMovie()
            filelist = [ f for f in os.listdir(self.templateMovDataDirectory) if  f.endswith(".png")]
            for f in filelist:
                os.remove(''.join([self.templateMovDataDirectory,'/',f]))
        
    def saveScreenShot(self):
        '''
        save screenshot to disk
        '''
        try: os.mkdir(self.screenShotDataDirectory)
        except: pass
        self.switch_to()
        pyglet.image.get_buffer_manager().get_color_buffer().save(''.join([self.screenShotDataDirectory,'/',self.networkName,'_',str(self.screenShotNumber),'.png']))
        self.screenShotNumber = self.screenShotNumber +1 
             
        
    def createMovie(self):
        '''
        open all movie-template-*.png's and create a movie out of it
        '''
        print "createMovie(): writing image data"
        
        frameSizeImage = read_png(''.join([self.templateMovDataDirectory,'/.screenShot',str(0),'.png']))
        frameSize = (np.shape(frameSizeImage)[1],np.shape(frameSizeImage)[0])
        
        try:   FFMpegWriter = animation.writers['mencoder']
        except: print "ERROR: Visualisation3D.createMovie(): mencoder libary is not installed, could not create movie!"; return
                
        try:
            fileName = ''.join([self.movieSaveDirectory,'/',self.networkName,'_',str(self.movieNumber),self.movieFileType])
            imageName = ''.join(['mf://',self.templateMovDataDirectory,'/.screenShot%d.png'])
            imageType = ''.join(['type=png:w=',str(frameSize[0]),':h=',str(frameSize[1]),':fps=24'])
            command = ('mencoder',
                       imageName,
                       '-mf',
                       imageType,
                       '-ovc',
                       'lavc',
                       '-lavcopts',
                       'vcodec=mpeg4',
                       '-oac',
                       'copy',
                       '-o',
                       fileName)
            os.spawnvp(os.P_WAIT, 'mencoder', command)
            self.movieNumber = self.movieNumber+1
            print "createMovie(): created movie sucessfull"
        except:
            print "ERROR: Visualisation3D.createMovie(): mencoder libary is not installed, could not create movie!"; return
               
        
        
    #--- functions for the visualisation speed
    def evaluateUpdateRate(self):
        '''
        evaluate the visualisation update rate and the timeStepIncrement
        to get realtime visualisation i.e. synchronise update times
        '''
        maxFrameRate = 60.
        minFrameRate = 24.
        numberOfUpdates = 50.
        
        dtSim     = self.vascularNetwork.dt
        timeSteps = len(self.vessels3D[self.vessels3D.keys()[0]].Psol)-1 # number of dt is one less then number of solutions
        totalTime = timeSteps*dtSim
        
        # approximate frame rate
        timeSolverSolve = []
        self.timeStepCurrent = int(self.timeStepsTotal/2.)
        for i in xrange(int(numberOfUpdates)):
            timeSolverSolveStart = time.clock()
            self.updateVisualisation(0.0)
            timeSolverSolve.append(time.clock()-timeSolverSolveStart)
         
        approximatedFrameRateUpdate = (1./(sum(timeSolverSolve[2::])/(numberOfUpdates-2)))-1
        
        if maxFrameRate < approximatedFrameRateUpdate: minFrameRate = maxFrameRate
        else: minFrameRate = approximatedFrameRateUpdate
                        
        x = 1./(dtSim*minFrameRate)
        x = int(np.ceil(x))+1
        frameRate = 1.0/(dtSim*x)
        
        self.numberOfCalls     = timeSteps/x
        self.timeStepIncrement = x
        self.timeStepIncrementRealTime = x
        self.updateTime = totalTime/self.numberOfCalls
        self.totalTime = totalTime
        
        # reset initial state
        self.timeStepCurrent = 0
        self.updateVisualisation(0.0)
        self.timeStepCurrent = self.timeStepIncrement    
    
    def visualisationSpeedHalfTime(self):
        self.timeStepIncrement = int(self.timeStepIncrementRealTime/2.)
    
    def visualisationSpeedRealTime(self):
        self.timeStepIncrement = int(self.timeStepIncrementRealTime)
        
    def visualisationSpeedUp(self):
        self.timeStepIncrement = min(int(self.timeStepIncrement+1),int(self.timeStepIncrementRealTime))
    
    def visualisationSpeedDown(self):
        self.timeStepIncrement = max(int(self.timeStepIncrement-1),1)
    
    def setVisualisationSpeed(self, value):
        self.timeStepIncrement = max(int(self.timeStepIncrementRealTime*value),1)
    
    #--- functions for the area movement
    def setAreaFactor(self,value):
        self.areaFactor[0] = np.min([np.max([1,value*self.areaFactorMaximum]),self.areaFactorMaximum])
        
    def enableWallMovement(self):
        self.wallMovement[0] = [True,False][self.wallMovement[0]]
    
    def findMaximumAreaFactor(self):
        '''
        Function to find the maximum area factor
        that r is allways greater then 0.2 * r0
        '''
        rfmax = 5000
        for vessel in self.vascularNetwork.vessels.itervalues():
            A0 = vessel.Asol[0]
            devisor = ((vessel.Asol[1::]-A0))
            devisor[devisor==0] = 1
            rf = ((0.1*vessel.Asol[1::]-A0))/ devisor
            ## take out values close to -inf
            rf[rf<0] = 5000
            rf = np.min(rf)           
            rfmax = np.min([rfmax,rf])
        self.areaFactorMaximum = rfmax
        self.areaFactor[0] = 1
                
    #--- functions for the look up tables  
    def changeColorTable(self):
        self.cLUT.changeColorTable()
    
    def changeQuantityLUT(self, quantity):
        self.quantityLUT[0] = quantity
        self.cLUT.changeLUTquantity(quantity)
    
    def enableWaveSplitLookUp(self):
        self.waveSplit[0] = [True,False][self.waveSplit[0]]
        self.cLUT.enableWaveSplitLookUp()
                        
    def enableLeftRightColoring(self):
        self.cLUT.enableLeftRightColoring()
    
    def findMinMaxPressure(self):
        maximaP = [1.e20, 0]
        maximaQ = [1.e20, 0]
        for vessel in self.vascularNetwork.vessels.itervalues():
            # maximum
            maximaP[1] = np.max([maximaP[1], np.max(vessel.Psol)])
            maximaQ[1] = np.max([maximaQ[1], np.max(vessel.Qsol)])
            # minimum
            maximaP[0] = np.min([maximaP[0], np.min(vessel.Psol)])
            maximaQ[0] = np.min([maximaQ[0], np.min(vessel.Qsol)])
            
        self.quantitiyMaxima = {'Pressure':maximaP, 'Flow':maximaQ}
            
         
class Vessel3D(Vessel):
    '''

    extended class of class Vessel() with functions for 3D manipulation and visualisation
    '''
    def __init__(self, vessel, batch, nPointsPerCircle, wallMovement, areaFactor, cLUT, quantityLUT, waveSplit):
        # init super class
        Vessel.__init__(self, vessel.Id, vessel.name)
        #pyglet.window.Screen.__init_(self,size)
        # apply all variables from super class instance to this class instance        
        
        self.oldVessel = vessel
        
        self.update(vessel.getVariableDict())   
        
        self.initialize({})
        self.quiet = True
        
       
        # new variables
        self.areaFactor     = areaFactor
        self.wallMovement   = wallMovement
        
        self.nPointsPerCircle = nPointsPerCircle
        
        self.cosSinZarray       = 0
        self.verticesInitial    = 0
        self.normalIndicesField = 0
        self.vertexList         = 0
        self.vertexList2        = 0
        self.batch              = batch
        
        self.viewNormals = False
        
        # LUT 
        self.cLUT = cLUT
        self.quantityLUT = quantityLUT
        self.waveSplit = waveSplit
        self.createWaveSplitSolutions()
        # functions
        self.createInitial3dVertice()
        self.waveSplitRange = {'Pressure': [[np.min(self.PsolF),np.max(self.PsolF)],
                                            [np.min(self.PsolB),np.max(self.PsolB)]],
                                'Flow':    [[np.min(self.QsolF),np.max(self.QsolF)],
                                            [np.min(self.QsolB),np.max(self.QsolB)]]}
            
    def createWaveSplitSolutions(self):
        '''
        calculate the forward and backward waves in advance
        '''
        numberOfTimeSteps = len(self.Psol)
        
        self.PsolF = np.zeros_like(self.Psol)
        self.PsolB = np.zeros_like(self.Psol)
        self.QsolF = np.zeros_like(self.Psol)
        self.QsolB = np.zeros_like(self.Psol)
        
        for n in xrange(int(self.N)):
            pf,pb,qf,qb =  linearWaveSplitting(self.Psol[:,[n]],self.Qsol[:,[n]],self.Asol[:,[n]],self.csol[:,[n]],self.rho)
            self.PsolF[1::,[n]] = pf.reshape(numberOfTimeSteps-1,1)
            self.PsolB[1::,[n]] = pb.reshape(numberOfTimeSteps-1,1)
            self.QsolF[1::,[n]] = qf.reshape(numberOfTimeSteps-1,1)
            self.QsolB[1::,[n]] = qb.reshape(numberOfTimeSteps-1,1)
        
    def createInitial3dVertice(self):            
        '''
        
        '''
        # creates verices
        nPointsPerCircle = self.nPointsPerCircle
        ## get properties of vessel
        N = int(self.N)
        z = self.z
        
        angleBetweenPoints  = 2*pi/nPointsPerCircle # needed that we do not have 2 points at start and end
        angleShift          = (pi+angleBetweenPoints)/2.0 # needed for left right coloring
        theta = np.linspace(0.0+angleShift, 2*pi-angleBetweenPoints+angleShift, nPointsPerCircle)
        ## set initial radius
        try:
            r = np.sqrt(self.Asol[0]/np.pi)  # try using solution data
        except: 
            r = np.sqrt(self.AsVector/np.pi) # use As vector if no solutionData is available 
            print "WARNING: 3dVisialisation.class3DVessel,createInitial3dVertice(): area solution available, taking As instead!"
                 
        ## create matrix with constants for cos,sin and z values
        cosRepeated = np.repeat([np.cos(theta)],N,axis=0).ravel()
        sinRepeated = np.repeat([np.sin(theta)],N,axis=0).ravel()
        zRepeated = np.repeat(z,nPointsPerCircle)
        self.cosSinZarray = np.array([cosRepeated,sinRepeated,zRepeated]).T
        ## create radius array and multiply it with self.cosSinZarray
        rrepeated = np.repeat(r,nPointsPerCircle).ravel()
        radiusArray = np.array([rrepeated,rrepeated,np.ones(nPointsPerCircle*N)]).T # needed for update
        vertices = self.cosSinZarray*radiusArray # needed for update
        ## rotate with rotation matrix to global system
        vertices = np.dot(vertices,self.rotToGlobalSys[0])
        ## move to proper position and ravel array
        vertices = (vertices+ self.positionStart[0].T)
        
        self.verticesInitial = np.copy(vertices)
        
        indeces = []        
        for i in xrange((N-1)*nPointsPerCircle):            
            if (i+1)%(nPointsPerCircle) == 0:
                indeces.extend([i,i-nPointsPerCircle+1,i+1,i+nPointsPerCircle])
            else:
                indeces.extend([i,i+1,i+nPointsPerCircle+1,i+nPointsPerCircle])
                
        ## create normals
        normals = np.empty_like(vertices)
        
        #self.normalIndicesOld = []
        normalIndicesStart1 = []
        normalIndicesStart2 = []
        normalIndicesField1 = []
        normalIndicesField2 = []
        normalIndicesField3 = []
        normalIndicesEnd1   = []
        normalIndicesEnd2   = []
        # create a normal for each vertex point
        numberOfVertices = len(vertices)
        for i in xrange(numberOfVertices):
            # find indices of vectors needed for the normals
            ip1 = i+1
            imP = i-nPointsPerCircle
            im1 = i-1
            ipP = i+nPointsPerCircle
            
            if (i)%(nPointsPerCircle) == 0:
                im1 = i+nPointsPerCircle-1
            
            if (i+1)%(nPointsPerCircle) == 0:
                ip1 = i-nPointsPerCircle+1              
                                    
            ## calculate normals
            if i < nPointsPerCircle: # at start: imP is not existing
                normals[i] = np.sum(np.cross(vertices[[ipP, ip1]]-vertices[i],vertices[[im1, ipP]]-vertices[i]),axis=0)
                #self.normalIndicesOld.append([i,[ipP, ip1],[im1, ipP]])
                
                normalIndicesStart1.extend([im1, ipP, ip1])
                normalIndicesStart2.extend([i,i,i])
                
            elif i >= (len(vertices)-nPointsPerCircle): #at end ring: ipP not existing
                normals[i] = np.sum(np.cross(vertices[[imP, im1]]-vertices[i],vertices[[ip1, imP]]-vertices[i]),axis=0)
                #self.normalIndicesOld.append([i,[imP, im1],[ip1, imP]])
                
                normalIndicesEnd1.extend([im1, ip1, imP])
                normalIndicesEnd2.extend([i,i,i])
                
            else: # in the middle
                normals[i] = np.sum(np.cross(vertices[[ipP, ip1, imP, im1]]-vertices[i],vertices[[im1, ipP, ip1, imP]]-vertices[i]),axis=0)
                #self.normalIndicesOld.append([i,[ipP, ip1, imP, im1],[im1, ipP, ip1, imP]])
                
                normalIndicesField1.extend([ipP, ip1, imP, im1])
                normalIndicesField2.extend([im1, ipP, ip1, imP]) 
                normalIndicesField3.extend([i,i,i,i])
          
        
        self.normalIndicesField = [normalIndicesField1,normalIndicesField2,normalIndicesField3]
        
          
        norm = np.linalg.norm(normals, axis= 1).reshape(len(normals),1)
        normals = normals /norm
        
        if self.viewNormals == True:
            ## create normal lines to visualize them
            normalEndpoint = vertices+(normals*0.2)
            normalLines = np.append(vertices.ravel(),normalEndpoint.ravel())
            numberOfNormals = len(vertices)
            normalLineIndeces = []
            # create normal indeces for the normal line
            for i in xrange(numberOfNormals):
                normalLineIndeces.extend([i,i+numberOfNormals])
        
        ## create color map
        # get color scalars depending on look up table for each grid node
        colorSolution = self.cLUT.getColors(self.Psol[0])
        # create RGB color vector for each vector of the vessel geometry
        color = np.repeat(colorSolution,nPointsPerCircle,axis=0).ravel()
                        
        vertices = vertices.ravel()
        normals  = normals.ravel()
                
        self.vertexList = self.batch.add_indexed(len(vertices)/3,  pyglet.gl.GL_QUADS, None, indeces,
                                                 ('v3f/stream',vertices),
                                                 ('n3f/stream',normals),
                                                 ('c3B/stream',color))
        
        if self.viewNormals == True:
            self.vertexList2 = self.batch.add_indexed(len(normalLines)/3,  pyglet.gl.GL_LINES, None,normalLineIndeces,
                                              ('v3f/stream',normalLines))
        
      
    def update3Dposition(self,timeStepCurrent, positionEndMother, rotToGlobalSysMother):
        '''
        
        '''
        nPointsPerCircle = self.nPointsPerCircle
        ## get properties of vessel
        N = int(self.N)
        ## 0.0 get updated rotational martices // or update them if movment is on
        if self.wallMovement[0]:
            ## 1.0 calculate new vertices
            areaFactor = self.areaFactor[0]
            r = np.sqrt( (self.Asol[0]*(1-areaFactor)+self.Asol[timeStepCurrent]*areaFactor) /np.pi)
        else: 
            areaFactor = self.areaFactor[0]
            r = np.sqrt(self.Asol[0]/np.pi)
        rrepeated = np.repeat(r,nPointsPerCircle).ravel()
        radiusArray = np.array([rrepeated,rrepeated,np.ones(nPointsPerCircle*N)]).T # needed for update
        vertices = self.cosSinZarray*radiusArray # needed for update
        ## rotate with rotation matrix to global system
        vertices = np.dot(vertices,self.rotToGlobalSys[timeStepCurrent])
        ## move to proper position and ravel array
        vertices = (vertices+ self.positionStart[timeStepCurrent].T)
        # update vertices
        self.vertexList.vertices = vertices.ravel()
                
        ## 2.0 calculate new normals
        ## create normals
        # calcualted normals of faces i.e. cross products
        crossed = np.cross(vertices[self.normalIndicesField[0]]-vertices[self.normalIndicesField[2]],vertices[self.normalIndicesField[1]]-vertices[self.normalIndicesField[2]])
        # sum normals
        normals = np.sum(crossed.reshape(nPointsPerCircle*(N-2),4,3),axis=1)
        # normalize
        norm = np.linalg.norm(normals, axis= 1).reshape(len(normals),1)        
        try:  normals = normals /norm
        except: pass
        # update normals
        self.vertexList.normals[nPointsPerCircle*3:nPointsPerCircle*(N-1)*3] = normals.ravel()
        
        
        leftRightColoringFactor = 1 # split coloring in left/right colorSplitFactor = 2 else colorSplitFactor = 1
        ## 3.0 calculate new colors
        ## create color map
        if self.quantityLUT[0] == 'Pressure':
        # get color scalars depending on look up table for each grid node
            if self.waveSplit[0] == False:
                colorSolution = self.cLUT.getColors(self.Psol[timeStepCurrent])
            else:
                colorSolution,leftRightColoringFactor = self.cLUT.getColorsWaveSplit(self.PsolF[timeStepCurrent],self.PsolB[timeStepCurrent],self.waveSplitRange)
        
        elif self.quantityLUT[0] == 'Flow':
        # get color scalars depending on look up table for each grid node
            if self.waveSplit[0] == False:
                colorSolution = self.cLUT.getColors(self.Qsol[timeStepCurrent])
            else:
                colorSolution,leftRightColoringFactor = self.cLUT.getColorsWaveSplit(self.QsolF[timeStepCurrent],self.QsolB[timeStepCurrent],self.waveSplitRange)
                
        
        # create RGB color vector for each vector of the vessel geometry
        color = np.repeat(colorSolution,nPointsPerCircle*leftRightColoringFactor,axis=0).ravel()
        # update color
        self.vertexList.colors = color        
                
        ## 5.0 update normal vertices
        if self.viewNormals:
            ## create normal lines to visualize them
            normalEndpoint = vertices[nPointsPerCircle:nPointsPerCircle*(N-1)]+(normals*0.2)
            normalLines = np.append(vertices[nPointsPerCircle:nPointsPerCircle*(N-1)].ravel(),normalEndpoint.ravel())
            self.vertexList2.vertices[nPointsPerCircle*3:nPointsPerCircle*(N-1)*3] = normalLines.ravel()
      

if __name__ == '__main__':
           
    optionsDict = parseOptions(['f','n','c'], visualisationOnly = True)
    
    networkName           = optionsDict['networkName']
    dataSetNumber         = optionsDict['dataSetNumber']
    connect               = optionsDict['connect']
             
    if networkName == None:
        #Visualisation3DGUI
        Visualisation3D()
        
    else:

        visualisation3D = Visualisation3D(networkName, dataSetNumber, connect)
            
    
        pyglet.app.run()  
        print "and exit()"
        pyglet.app.exit()
        ### garbage collection
        gc.collect()
        del visualisation3D
    