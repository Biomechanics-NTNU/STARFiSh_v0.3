import pprint
from copy import deepcopy

from moduleGrids import *
from classCompliance import *

import sys, os
cur = os.path.dirname( os.path.realpath( __file__ ) )
sys.path.append(cur+'/../'+'/UtilityLib')
from constants import vesselComplianceElements

from math import pi, cos, sin
import numpy as np
from copy import deepcopy

__author__ = "Vinzenz Gregor Eck"
__version__ = "0.3"

class Vessel(object):
    '''
    Class representing a vessel in a vascular Network
    '''
    number = 0 
    quiet = False
    
    def __init__(self, Id= None, name = None):
        '''
        Constructor
        '''
        ### Input properties
        
        ## attributes
        self.Id     = Id                              # id of the vessel
        self.name   = name                            # name of the vessel
        
        ## topology properties
        self.startNode          = None              # node id of starting node  
        self.endNode            = None              # node id of ending node            
        self.leftDaughter       = None              # id of left daugther vessel
        self.rightDaughter      = None              # id of right daugther vessel
        self.leftMother         = None              # id of left mother
        self.rightMother        = None              # id of right mother
        self.angleXMother       = 0                 # rotation angle in rad around x-axis relative to mother vessel
        self.angleYMother       = 0                 # rotation angle in rad around y-axis relative to mother vessel
        self.angleZMother       = 0                 # rotation angle in rad around z-axis relative to mother vessel
                
        ## geometry properties
        self.geometryType       = 'uniform'         # geommetry type: 'Uni', 'Cone', 'Cons'
        self.radiusProximal     = 0.01              # radius at vessel begin // in hole Vessel if radiusB = 0
        self.radiusDistal       = 0.01              # radius at vessel end
        self.length             = 0.1               # length of the vessel
        self.N                  = 5                 # number of gridpoints
        
        ## compliance properties
        self.complianceType     = 'Hayashi'         # type of the compliance see classCompliance
        self.constantCompliance = False             # use constant compliance C0 calculated with Ps for hole simulation
        self.externalPressure   = 0.                # external pressure of the vessel usually 0 for head arteries 2000 Pa
        self.Ps                 = 0.                # pressure at state S
        self.As                 = None              # Area state S input from XML if None it is calculate from given radius
        
        self.wallThickness      = None              # wallthickness
        self.youngModulus       = None              # youngModulus
        self.distensibility     = None              # arterial distensibility for Reymond compliance model
        self.betaExponential    = None              # material Parameter of Exponential compliance model
        self.betaHayashi        = 1.                # material Parameter of Hayashi compliance model
        self.betaLaplace        = None              # material Parameter of Laplace compliance model
        
        # FLUID properties 
        self.applyGlobalFluid   = True              # bool: apply global fluid properties or the ones stored in vessel XML
        self.my                 = 1.e-6             # blood viscosity
        self.rho                = 1050.             # blood density
        self.gamma              = None              # velocity profile gamma correction                              
        
        # vascular Polynomial chaos data dictionary with uncertainty variables
        self.polyChaos          = {}                # dict with uncertain variables { 'variableName': [dist coeff coeff]}
        
        ### Calculated properties
        # positions in space 
        # realtive to global system
        self.positionStart  = np.zeros((1,3))         # instanteanoues position of vessel start point in the global system
        self.positionEnd    = np.zeros((1,3))         # instanteanoues position of vessel start point in the global system
        self.rotToGlobalSys = np.array([np.eye(3)])   # rotational matrix to the global system
        
        ## motion
        self.angleXMotherTime   = []                # describes angular motion over time; set up by flow solver
        self.angleYMotherTime   = []                # describes angular motion over time; set up by flow solver
        self.angleZMotherTime   = []                # describes angular motion over time; set up by flow solver
        
        # gravity 
        self.netGravity         = [0.0]             # net gravity on the vessel
        self.gravityConstant    = -9.81             # gravity constant of the gravitational field
        
        # GRID properties
        self.z                  = None              # axial coordinates of the vessel
        self.dz                 = None              # axial grid spacing
        self.A0                 = None              # initial Area
        self.AsVector           = None              # vector for the total grid containing area at state S 
        
        # SOLID properties   
        self.compliance         = None              # the compliance class reference
        self.P_Cexp             = None              # pressuer P0 for the ComplianceType exponential
        self.C                  = None              # Compliance function C(P) for all nodes
        self.C_nID              = None              # Compliance function C(P,NodeId) for specific Node with NodeID
        self.A                  = None              # Area function A(P) for all nodes
        self.A_nID              = None              # Area function A(P,NodeId) for specific Node with NodeID
        self.c                  = self.waveSpeed    # wave speed function c(rho,A,C) for all nodes /single if given single values
        
        self.resistance         = None              # vessel resistance
        self.womersleyNumber    = None              # Womersley number
        
        
        ## Solution Data
        self.Psol      = None  # pressure
        self.Qsol      = None  # flow 
        self.Asol      = None  # area 
        # calculated values for post processing
        self.csol      = None # wave speed
        self.vsol      = None # mean velocity
        
        self.quiet = False
        Vessel.number += 1
        
        
        ### old variables of xml version 3.0 very outdated to be deleted soon
#         self.geom = None                                        # geommetry type: 'Uni', 'Cone', 'Cons'
#         self.lengthInterval = None
#         self.radiusA = None                                  
#         self.radiusAInterval = None
#         self.radiusB = None                                     # radius at vessel end
#         self.radiusBInterval = None
#         self.NInterval = None                              
#         
#         self.Pfunc = None                                      # pressure function: 'True', 'False
#         self.wallThicknessInterval = None
#         self.youngModulusInterval= None
#         self.betaInterval = None
#         
#         self.betaVector = None
#         
#         self.beta               = None              # beta value from XML-file
#         
#         self.myInterval = None
#         self.rhoInterval = None 
                
    def initialize(self, globalFluid):
        '''
        Initialisation of the vessel.
        This method calculates and set up, the compliance, gird, fluid and resistance of the Vessel.
        '''
        # initialze fluid
        for key,value in globalFluid.iteritems():
            
            try:
                if self.applyGlobalFluid == True:
                    self.__getattribute__(key)
                    self.__setattr__(key,value)
            except:
                if key != 'venousPressure':
                    print "ERROR: classVessel.initialize(): Fluid initialisation could not update varibale {} of vessel {}!".format(key,self.Id)
        
        # initialize gird
        try:
            self.z,self.dz,self.A0 =  eval(self.geometryType)(self.length, self.radiusProximal, self.radiusDistal, self.N)
        except:
            print "ERROR: classVessel.initialize(): Grid calculation of vessel {}!".format(self.Id)
            exit()
        
        ## calculate compliance reference area vector As
        try:
            if self.geometryType == 'uniform':
                if self.As == None:
                    AsProximal = self.radiusProximal**2.0*np.pi
                    As_distal = AsProximal
                    self.AsVector = np.linspace(AsProximal,As_distal,self.N)
                    #print "WARNING: no reference Area for the Compliance As is given for vessel {} ! \n         As is now calculatet with given radiusA".format(self.Id)
                else:
                    self.AsVector = np.ones(self.N)*self.As
            elif self.geometryType == 'cone':
                AsProximal = self.radiusProximal**2.0*np.pi
                As_distal = self.radiusDistal**2.0*np.pi
                #AsAveraged = (AsProximal+As_distal)/2.0
                self.AsVector = np.linspace(AsProximal,As_distal,self.N)
                #self.AsVector = np.linspace(AsAveraged,AsAveraged,self.N)
                #print "WARNING: no reference Area for the Compliance As is given for vessel {} ! \n         As is now calculatet with given radiusA and radiusB".format(self.Id)      
            else:
                print "WARNING: classVessel.initialize(): no grid geometry not proper defined for vessel {} ! ".format(self.Id)      
        except:
            print "ERROR: classVessel.initialize(): As calculation of vessel {}!".format(self.Id)
            
        ## initialize compliance
        try:
            self.compliance = eval(self.complianceType)(self.rho,self.AsVector)
            # initialize compliance element
            complianceDict = {}
            for variable in vesselComplianceElements[self.complianceType]:
                if variable not in ['As','constantCompliance']:
                    complianceDict[variable] = self.getVariableValue(variable)
            self.compliance.initialize(complianceDict)
            
            # apply functionals
            self.A  = self.compliance.A
            self.A_nID = self.compliance.A_Node
            
            if self.constantCompliance == False:
                self.C     = self.compliance.C
                self.C_nID = self.compliance.C_Node
            else:
                self.C     = self.compliance.C0
                self.C_nID = self.compliance.C0_Node
        
        except:
            print "ERROR: classVessel.initialize(): init Compliance of vessel {}!".format(self.Id)
            exit()
                
        # set up Resistance
        try:
            # calculate peuisell resistance R = (8*mu*L) / (pi*r**4)
            if self.geometryType == "uniform":
                self.resistance = 8*self.my*self.length / (pi*self.radiusProximal**4)
            elif self.geometryType == "cone":
                self.resistance = 8*self.my*self.length / (pi*((self.radiusProximal+self.radiusDistal)/2)**4)

        except:
            print "ERROR: classVessel.initialize(): in calculating resistance of vessel {}!".format(self.Id)
        
        # calculate Womersley number
        
        
        #set hooks waveSpeed function
        if self.c == None: self.c = self.waveSpeed  
        
     
    def waveSpeed(self,Area,Compliance,sqrt= np.sqrt):
        '''
        Input:
            Area: Area array or float at one ID
            Compliance: Compliance array or float at one ID 
            NOTE: it does not check if len(Area) == len(Compliance)
        Output:
            waveSpeed : sqrt(A/(rho*C))
        '''
        try : return sqrt(Area/(self.rho*Compliance))
        except:
            print "ERROR: classVessel.waveSpeed(): wave speed calculation of vessel id {}!".format(self.Id)
            if Area.any() < 0: print "     Area < 0 !!!"
            if Compliance.any() < 0: print "     Compliance < 0 !!!"
            exit()
    
    def postProcessing(self):
        '''
        Post processing methode
        calulate wave speed vector from solution data for pressure and flow
        calcuate mean flow from the solution data of flow and area
        '''
        self.csol = self.waveSpeed(self.Asol,self.C(self.Psol))
        self.vsol = self.Qsol/self.Asol
        
    def Impedance(self,Pressure):
        '''
        Input:
            Pressure array for the hole grid
        Output:
            Impedance : 1.0/(c*C))
        '''
        Compliance = self.C(Pressure)
        Area = self.A(Pressure)
        c = self.waveSpeed(Area, Compliance)
        return 1.0/(c*Compliance) 
    
    def update(self, Dict):
        '''
        updates the vessel data using a dictionary in from of 
	    dataDict = {'variableName': value}
        '''
        for key,value in Dict.iteritems():
            try:
                self.__getattribute__(key)
                self.__setattr__(key,value)
            except:
                print "ERROR Vessel.updateData (vesselId {}) Wrong key: {}, could not update varibale".format(self.Id, key)
    
    def getVariableValue(self,variableName):
        '''
        Returns value of variable with name : variableName
        States Error if not such variable
        '''
        try:
            return self.__getattribute__(variableName)
        except: 
            print "ERROR Vessel.getVariable() : vessel has no variable {}".format(variableName)
                       
    def getVariableDict(self):
        '''
        Returns a deep copy of the class variable dict
        '''
        return deepcopy(self.__dict__)
    
    def printToConsole(self):
        '''
        writes all variables and their values to the console
        '''
        print "----------------"
        print "    Vessel %d"%self.Id,"\n"
        for variable,value in self.__dict__.iteritems():
            try: 
                print " {:<20} {:>8}".format(variable,value)
            except: print " {:<20} {:>8}".format(variable,'None')
            
    def caculatePositionAndGravity(self, n, positionEndMother, rotToGlobalSysMother):
        '''
        calculates the position and gravity for given time point n and inserts these
        in the corresponding vectos
                
        Initializing net gravity on the vessels.
        '''
         
        try:    angleXMother = self.angleXMotherTime[n]
        except: angleXMother = self.angleXMother
        
        try:    angleYMother = self.angleYMotherTime[n]
        except: angleYMother = self.angleYMother
        
        try:    angleZMother = self.angleZMotherTime[n]
        except: angleZMother = self.angleZMother
        
        rotToGlobalSys = rotToGlobalSysMother
        
        # 2. assamble rot matrices,
        # 2.1 calculate local rotmatrix to mother depending on 
        ## x axis, angleXMother
        if angleXMother != 0:
            rotDaughterMotherX = np.array([[1, 0                     , 0                     ],
                                           [0, cos(angleXMother), sin(angleXMother)],
                                           [0,-sin(angleXMother), cos(angleXMother)]])
            
            rotToGlobalSys = np.dot(rotDaughterMotherX,rotToGlobalSys)    
            
        ## y axis, angleYMother
        if angleYMother != 0:
            rotDaughterMotherY = np.array ([[ cos(angleYMother), 0,  -sin(angleYMother)],
                                            [ 0                     , 1,  0                      ],
                                            [ sin(angleYMother), 0,  cos(angleYMother) ]])
            
            rotToGlobalSys = np.dot(rotDaughterMotherY,rotToGlobalSys)
            
        ## z axis, angleZMother
        if angleZMother != 0:
            rotDaughterMotherZ = np.array ([[  cos(angleZMother), sin(angleZMother),  0],
                                            [ -sin(angleZMother), cos(angleZMother),  0],
                                            [  0,                      0,                       1]])
            
            rotToGlobalSys = np.dot(rotDaughterMotherZ,rotToGlobalSys)
            
        # 3. calulate pos end
        positionEnd = np.array([0,0,self.length])
        positionEnd = (np.dot(positionEnd,rotToGlobalSys) + positionEndMother)
        
        gravityVector = np.array([0,0,self.gravityConstant])
        netGravity = np.dot(gravityVector,rotToGlobalSys)[2] 
        
        self.positionStart[n]  = positionEndMother ## positionStart
        self.positionEnd[n]    = positionEnd
        self.rotToGlobalSys[n] = rotToGlobalSys
        self.netGravity[n]     = netGravity
        
    


