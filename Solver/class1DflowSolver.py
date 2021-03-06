import numpy as np
import scipy as sc
from itertools import islice

import sys,os
from classBoundaryConditions import VaryingElastance, Valve
from classVascularNetwork import VascularNetwork

# set the path relative to THIS file not the executing file!
cur = os.path.dirname( os.path.realpath( __file__ ) )
sys.path.append(cur+'/NetworkLib')

from classBoundarys import Boundary

from classSystemEquations import *
from classConnections import *
from classFields import *
from classCommunicators import *

import pprint 
from copy import copy as copy 

import gc

__author__ = "Vinzenz Gregor Eck"
__version__ = "0.3"

class FlowSolver(object):

    
    def __init__(self,vascularNetwork, quiet=False):
        '''
        Constructor       
        '''
        if vascularNetwork == None: print "ERROR: No vascularNetwork given!" / exit()
        assert isinstance(vascularNetwork, VascularNetwork)
        # the vascular network to solve
        self.vascularNetwork = vascularNetwork
        self.vascularNetwork.quiet = quiet
        
        self.vessels = self.vascularNetwork.vessels
        self.fields = {}
        # the boundarys of the network { vesselID : [<instance>::classBoundary_02(Characteristics.py), .. ]}
        # 1 boundary for each start/end-vessel except if only 1 vessel in the network
        self.boundarys = {}
        # the system Equations of each vessel { vesselID : <instance>::classSystemEquations(SystemEquations) }
        self.systemEquations = {}
        # the connections of the Network { motherVesselID : <instance>::classConnections, ...}
        self.connections  = {}                
        # the communicator objects of the vascularNetwork {communicatorID : <instance>::classCommunicator}
        self.communicators = {}
        # list of numerical objects (field,connection,boundary objects as in the traversing list)
        self.numericalObjects = []
        # time step
        self.dt = None
        # number of Timesteps
        self.Tsteps = None
        # total Simulation time
        self.totalTime = None
        
        # timestep Counter
        self.n = [0]
        
        # Set div output 
        self.output = {}
        
        self.simplifyEigenvalues = self.vascularNetwork.simplifyEigenvalues
        
        self.riemannInvariantUnitBase = self.vascularNetwork.riemannInvariantUnitBase
        
        self.solvingSchemeConnections = self.vascularNetwork.solvingSchemeConnections
       
        # bool for cfl meshing with automatic grid adaption
        self.automaticGridAdaptation = self.vascularNetwork.automaticGridAdaptation
        
        # rigidAreas True: A = A0 False: A = A(P)
        self.rigidAreas = self.vascularNetwork.rigidAreas
        
#         # Define the output of A, dependend on the characteristic system 0.1
#         if self.rigidArea == '02':
#             self.AFunction = self.AFunctionSys0
#         elif self.rigidArea == '2': 
#             self.AFunction = self.AFunctionSys1

#         self.solvingScheme = self.vascularNetwork.networkSolver['solvingScheme'] 
#         if self.solvingScheme == 'MacCormack':
#             self.solve = self.MacCormack  
#         elif self.solvingScheme == 'MacCormack_Field': 
#             self.solve = self.MacCormack_Field
        
        #define solve function
        self.solve = self.MacCormack_Field
        
        ## for future use !!
        self.multiProcessing = False
                
        # initialize system
        self.vascularNetwork.initialize(initializeForSimulation = True)
        
        self.initializeTimeVariables(quiet)
        
        self.initializeSolutionMatrices()
        self.initializeSystemEquations()
        self.initializeBoundarys()
        self.initializeConnections()
        self.initializeFields()
        self.initializeCommunicators()
        self.initializeNumericalObjectList()
        if quiet==False:
            self.initOutput() # feedback
        self.quiet = quiet 
       
        

    '''       
    ########################################################################################
    # initialisation Methods
    ########################################################################################
    '''
        
    
    def calcTimeStep(self,dz,c,CFL):
        return (CFL*dz)/c

    def initializeTimeVariables(self, quiet):
        '''
        initialize time variable dt and Tstep
        '''
        self.totalTime = self.vascularNetwork.totalTime
                
        initialValues = self.vascularNetwork.initialValues
        
        dt_min,dz_min,c_max,gridNodens = [],[],[],[]
        #create waveSpeed Log file
        logfileData = {}
        
        dt = self.totalTime
        for vessel in self.vessels.itervalues():
        # Calculate time variables
            #estimate initial pressure
            p0,p1 = initialValues[vessel.Id]['Pressure']
            
            initialPressure = np.linspace(p0,p1,vessel.N) 
            
            A0_max = max(vessel.A(initialPressure))
            #c_high = vessel.c(A0_max,vessel.initialPressure)
            
            #c_high1 = vessel.c(A0_max,0)
            Compliance = vessel.C(initialPressure)
            c_high = vessel.c(A0_max, Compliance)
            
            dz_low = min(vessel.dz)
            dt = self.calcTimeStep(dz_low,c_high,self.vascularNetwork.CFL) 
            c_max = np.append(c_max,c_high)
            dt_min = np.append(dt_min,dt)
            dz_min = np.append(dz_min,dz_low)
            gridNodens = np.append(gridNodens,vessel.N)
            
            logfileData[vessel.Id] = [max(c_high),min(c_high),min(dt),vessel.dz,vessel.N]
            
        # Set time variables 
        self.dt = min(dt_min)

        self.Tsteps = int(round(self.totalTime/self.dt))

        self.output['dz_min']     = min(dz_min)
        self.output['c_min']      = min(c_max)
        self.output['c_max']      = max(c_max)
        self.output['gridNodens'] = sum(gridNodens)
        
        self.output['CFLcorrect'] = []
        
        automaticGridCorrection = {}
        
        logfile = open(str(cur+'/../'+'LOGcurrentWaveSpeed.txt'),'wb')
        logfile2 = open(str(cur+'/../'+'LOGproposedGrid.txt'),'wb')
        CFL = self.vascularNetwork.CFL
        for vesselT,data in logfileData.iteritems():
            #number of deltaX
            Nnew = int((sum(data[3])*CFL/(self.dt*data[0])))
            
            L = sum(data[3])
            #calculate number of dx: M = L/dx
            Mnew = int(L*CFL/self.dt/data[0])
                        
            dz_new = L/Mnew
            dt_new = self.calcTimeStep(dz_new,data[0],CFL)
            
            while dt_new < self.dt:
                Mnew = Mnew-1
                dz_new = L/Mnew
                dt_new = self.calcTimeStep(dz_new,data[0],CFL)
            
            # calculate gridpoints N = N+1
            Nnew = Mnew+1
            if Nnew == data[4]: 
                logfile.write(''.join(['vessel ',str(vesselT), ' c_max: %2.3f'%(data[0]),'   dt (ms) %2.6f'%(data[2]*1.0E3),' | res CFL: %2.3f'%(data[0]*self.dt/min(data[3])),' || already best N', '\n']))
                #self.output['CFLcorrect'].append(' '.join([str(vesselT).rjust(2).ljust(15),' no correction']))
            else: 
                logfile.write(''.join(['vessel ',str(vesselT), ' c_max: %2.3f'%(data[0]),'   dt (ms) %2.6f'%(data[2]*1.0E3),' | res CFL: %2.3f'%(data[0]*self.dt/min(data[3])),' || prop: N%2.0f'%Nnew,'    dN %.0f'%(Nnew-data[4]),'    dtNew %.4f'%(dt_new*1.e3),'   CFL: %2.3f'%(data[0]*self.dt/dz_new), '\n']))
                self.output['CFLcorrect'].append(' '.join([str(vesselT).ljust(3),'|',str(int(data[4])).ljust(3),'->',str(Nnew).ljust(3),'|', '%2.3f'%(data[0]*self.dt/min(data[3])),'->', '%2.3f'%(data[0]*self.dt/dz_new),'|']))
                automaticGridCorrection[vesselT] = Nnew
            logfile2.write(''.join([str(int(Nnew)),'\n']))
        logfile.close()
        logfile2.close()    
        self.vascularNetwork.dt = self.dt
        
        if self.output['CFLcorrect'] != []:
            if quiet == False:
                print '====================================='
                print '___CFL-correction: grid adaptation___'
                print 'Id  | gridPoints |       CFL      |'
                print '    | now -> new |   now -> new   |'
                for CFLcorr in self.output['CFLcorrect']:
                    print '%s' % (CFLcorr)
                
            if automaticGridCorrection != {}:
                gridCorrection = 'ohYesDoItPlease'
                if self.automaticGridAdaptation == True: gridCorrection = 'y'
                while  gridCorrection not in (' ','','y','Y','n','N'): 
                    gridCorrection = raw_input('Do you whish to adapt grid? (yes [<ENTER>,<y>,<Y>]/no [<n>,<N>])')
                
                if gridCorrection in (' ','','y','Y'):
                    #if quiet == False: print ' proceed with: grid aptation for vessels {} \n'.format(automaticGridCorrection.keys())
                    arrayN = 0
                    for vesselId,Nnew in automaticGridCorrection.iteritems():
                        self.vessels[vesselId].N = Nnew
                        self.vessels[vesselId].initialize({})
                        
                        newOutput = self.output['CFLcorrect'][arrayN].split(' no')[0]
                        self.output['CFLcorrect'][arrayN] = ' '.join([newOutput,'yes'])
                        arrayN = arrayN+1
        gridNodens = 0
        for vessel in self.vascularNetwork.vessels.itervalues():
            gridNodens += vessel.N 
            
        self.output['gridNodens'] = int(gridNodens)
        
        
    def initializeSolutionMatrices(self):
        '''
        initialize solution matrices
        '''
        
        initialValues = self.vascularNetwork.initialValues
        
        for vesselId,vessel in self.vessels.iteritems():
            
            vessel.Psol = np.ones((self.Tsteps,vessel.N))
            vessel.Qsol = np.zeros((self.Tsteps,vessel.N))
            vessel.Asol = np.zeros((self.Tsteps,vessel.N))
            
            
            try:
                p0,p1 = initialValues[vesselId]['Pressure']
                Qm    = initialValues[vesselId]['Flow']
                vessel.Psol[0] = np.linspace(p0,p1,vessel.N)   
                vessel.Qsol[0] = np.ones((1,vessel.N))*Qm  
            except:
                print "Error: cFS could not use initial values from network"
                pass
            vessel.Asol[0] = np.ones((1,vessel.N))*vessel.A(self.vessels[vesselId].Psol[0])   
                   
            vessel.positionStart    = np.zeros((self.Tsteps,3))
            vessel.positionEnd      = np.zeros((self.Tsteps,3))
            vessel.rotToGlobalSys   = np.zeros((self.Tsteps,3,3))
            vessel.netGravity       = np.zeros((self.Tsteps,1))
                       
        ## initialse varying elastance model
        for vesselId,boundaryConditions in self.vascularNetwork.boundaryConditions.iteritems():
            for bC in boundaryConditions:
                if isinstance(bC, VaryingElastance):
                    Qm    = initialValues[vesselId]['Flow']
                    bC.update({'aorticFlowPreviousTimestep':Qm})
                    bC.initializeSolutionVectors(self.Tsteps)
                    
        motionDict = {}         
        for vesselId,angleDict in motionDict.iteritems():
            self.vessels[vesselId].update(angleDict)
            
        ## calculate gravity and positions
        for n in xrange(self.Tsteps-1):            
            self.vascularNetwork.calculate3DpositionsAndGravity(n)
            
        ## calculate venous pressure for windkessel
        self.vascularNetwork.initializeVenousGravityPressureTime(self.Tsteps)
         
    def initializeSystemEquations(self):
        '''
        initialize system Equations
        '''
        for vesselId,vessel in self.vessels.iteritems():
            self.systemEquations[vesselId] = System(vessel,
                                                    self.simplifyEigenvalues,
                                                    self.riemannInvariantUnitBase,
                                                    self.n,
                                                    self.dt)
            # initialize system equations
            self.systemEquations[vesselId].updateSystem(self.vessels[vesselId].Psol[0],
                                                        self.vessels[vesselId].Qsol[0],
                                                        self.vessels[vesselId].Asol[0])
         
    def initializeBoundarys(self):
        '''
        initialize boundarys
        '''
        if len(self.vessels) == 1:
            rootId = self.vascularNetwork.root
            bcList0 = []
            bcList1 = []
            for bc in self.vascularNetwork.boundaryConditions[rootId]:
                if bc.position == 0:     bcList0.append(bc)
                elif bc.position == -1:  bcList1.append(bc)
            self.boundarys[rootId] = [Boundary( self.vessels[rootId],
                                                bcList0,
                                                self.rigidAreas,
                                                self.dt,
                                                self.n,
                                                self.Tsteps,
                                                self.systemEquations[rootId]),
                                      Boundary( self.vessels[rootId],
                                                bcList1,
                                                self.rigidAreas,
                                                self.dt,
                                                self.n,
                                                self.Tsteps,
                                                self.systemEquations[rootId])]
            self.output['BndrNR'] = 2
        else:
            for vesselId,boundaryConditions in self.vascularNetwork.boundaryConditions.iteritems():
                self.boundarys[vesselId] = [  Boundary( self.vessels[vesselId],
                                                        boundaryConditions,
                                                        self.rigidAreas,
                                                        self.dt,
                                                        self.n,
                                                        self.Tsteps,
                                                        self.systemEquations[vesselId])]
                
            self.output['BndrNR'] = len(self.boundarys)
    
    def initializeConnections(self):
        '''
        initialize Connections of the network
        by traversing the network tree
        '''
        treeList = self.vascularNetwork.treeTraverseList
        
        for leftMother,rightMother,leftDaughter,rightDaughter  in self.vascularNetwork.treeTraverseConnections:  
            ## link
            if rightMother == None and rightDaughter == None:
                self.connections[leftMother] = Link(  self.vessels[leftMother],
                                                      self.systemEquations[leftMother],
                                                      self.vessels[leftDaughter],
                                                      self.systemEquations[leftDaughter],
                                                      self.n,
                                                      self.dt,
                                                      self.rigidAreas,
                                                      self.solvingSchemeConnections)
            ## bifurcation
            elif rightMother == None:
                self.connections[leftMother] = Bifurcation(  self.vessels[leftMother],
                                                             self.systemEquations[leftMother],
                                                             self.vessels[leftDaughter],
                                                             self.systemEquations[leftDaughter],
                                                             self.vessels[rightDaughter],
                                                             self.systemEquations[rightDaughter],
                                                             self.n,
                                                             self.dt,
                                                             self.rigidAreas,
                                                             self.solvingSchemeConnections)
            ## anastomosis
            elif rightDaughter == None:
                anastomosisId = leftMother
                if treeList.index(leftMother) > treeList.index(rightMother):
                    anastomosisId = rightMother
                self.connections[anastomosisId] = Anastomosis(self.vessels[leftMother],
                                                             self.systemEquations[leftMother],
                                                             self.vessels[rightMother],
                                                             self.systemEquations[rightMother],
                                                             self.vessels[leftDaughter],
                                                             self.systemEquations[leftDaughter],
                                                             self.n,
                                                             self.dt,
                                                             self.rigidAreas,
                                                             self.solvingSchemeConnections)
        
    def initializeFields(self):
        
        for vesselId,vessel in self.vessels.iteritems():    
            self.fields[vesselId] = Field(  vessel,
                                            self.n,
                                            self.dt, 
                                            self.systemEquations[vesselId],
                                            self.rigidAreas,)
    
    def initializeCommunicators(self):
        
        
        #print 'cFS435 Communicators',self.vascularNetwork.communicators
        for comId, comData in self.vascularNetwork.communicators.iteritems():
                      
            ## for baro receptor and visualisation
            try:
                data = {'Pressure': self.vessels[comData['vesselId']].Psol,
                        'Flow'    : self.vessels[comData['vesselId']].Qsol,
                        'Area'    : self.vessels[comData['vesselId']].Asol
                        }
                comData['data']           = data
            except: pass
                                    
            comData['n']              = self.n
            comData['dt']             = self.dt
            
            self.communicators[comId] = eval(comData['comType'])(comData)
       
       
    def initializeNumericalObjectList(self):
        '''
        ## fill numObjectList (self.numericalObjects) traversing the treeList 
        # 1. add root boundary
        # 2  add vessels
        # 3  add connection or distal boundary condition
        # 4. repeat 2,3 for the hole tree 
        # 5. add communicators
        # 6. add blocking Wait if multiprocessing
        '''
        
        # get treetraversing list
        treeList = self.vascularNetwork.treeTraverseList
        
        singleVesselNetwork = False
        if len(self.vessels) == 1:
            singleVesselNetwork = True
        
        for vesselId in treeList:
            int(vesselId)
            ## check if root add BC
            try:
                if vesselId == self.vascularNetwork.root:
                    self.numericalObjects.append(self.boundarys[vesselId][0])
            except: pass
            
            ## add field
            self.numericalObjects.append(self.fields[vesselId])
            
            ## try add Connection
            try: self.numericalObjects.append(self.connections[vesselId])    
            except: pass
            
            ## try add distal BC
            try:
                if vesselId in self.vascularNetwork.boundaryVessels:
                    if singleVesselNetwork == False:
                        self.numericalObjects.append(self.boundarys[vesselId][0])
                    else:
                        self.numericalObjects.append(self.boundarys[vesselId][1])
            except: pass
        
        for communicator in self.communicators.itervalues():
            self.numericalObjects.append(communicator) 
            try:    communicator.startRealTimeVisualisation()
            except: pass
            
        if self.multiProcessing == True:
            self.numericalObjects.append("HOlD")
               
    def initOutput(self):
        '''
        initialize solution matrices
        '''
        #print '====================================='
        print '___________Time variables ___________'
        print '%-20s %2.1f' % ('totaltime (sec)',self.totalTime)
        print '%-20s %2.3f' % ('dt (ms)',self.dt*1.0E3)
        print '%-20s %4d' % ('Tsteps',self.Tsteps)
        print '___________Div variables ____________'
        print '%-20s %2.1f' % ('Q init (ml s-1)',self.vascularNetwork.initialValues[self.vascularNetwork.root]['Flow']*1.e6)
        print '%-20s %2.1f' % ('P init (mmHg)',self.vascularNetwork.initialValues[self.vascularNetwork.root]['Pressure'][0]/133.32)
        try: print '%-20s %2.1f' % ('R_cum (mmHg s ml-1)',self.vascularNetwork.Rcum[self.vascularNetwork.root]/133.32*1.e-6)
        except: pass
        print '%-20s %2.1f' % ('CFL init max',self.vascularNetwork.CFL)
        print '%-20s %2.1f' % ('dz min (mm)',self.output['dz_min']*1.0E3)
        print '%-20s %2.1f' % ('c min (m/s)',self.output['c_min'])
        print '%-20s %2.1f' % ('c max (m/s)',self.output['c_max'])
        print '%-20s %4d' % ('Grid nodens',self.output['gridNodens'])
        print '___________Num variables ____________'
        print '%-20s %4d' % ('NumObj',len(self.fields)+len(self.boundarys)+len(self.connections))
        print '%-20s %4d' % ('NumFields',len(self.fields))
        print '%-20s %4d' % ('NumConnections',len(self.connections))
        print '%-20s %4d' % ('NumBoundarys',len(self.boundarys))
        print '%-20s %4d' % ('NumCommunicators',len(self.communicators))
        print '%-20s %4d' % ('NumObj calls',len(self.numericalObjects)*self.Tsteps)                 
        print '===================================== \n'
        
            
    '''
    ########################################################################################
    # Solver Methods:
    #
    #    MacCormack_Field
    #    
    #
    ########################################################################################
    '''
    
    def MacCormack_Field(self):
        '''
        MacCormack solver method with forward-euler time steping,
        Using either Characteristic system 0 or 1 as defined in the XML-file.
        
        This method is solving the system by looping through the defined network
        imposing the boundary conditions based on Riemann Invariants and then solving the vessels, 
        conncetions, bifucations with a predictor-corrector step method
        '''
        print "Solving system ..."
                
        # original
        for n in xrange(self.Tsteps-1):
            
            self.n[0] = n
            for numericalObject in self.numericalObjects:
                numericalObject()
                                                    
            
        print "\nSystem solved!"
        
        if self.quiet == False:
            
            print '\n====================================='
            print '___Blood volume - consistency check___'
            
            print 'Vessels  init (ml)  sol (ml)  diff'
            vesselsInit = 0
            vesselsSol  = 0
            for vesselId,vessel in self.vessels.iteritems():
                A1 = self.vessels[vesselId].Asol[0][0:-1]
                A2 = self.vessels[vesselId].Asol[0][1:]
                volumeInit = np.sum(vessel.dz*(A1+A2+np.sqrt(A1*A2))/3.0)*1.e6
                A1 = self.vessels[vesselId].Asol[-1][0:-1]
                A2 = self.vessels[vesselId].Asol[-1][1:]
                volumeSol  = np.sum(vessel.dz*(A1+A2+np.sqrt(A1*A2))/3.0)*1.e6
                diff = volumeInit-volumeSol
                vesselsInit += volumeInit
                vesselsSol  += volumeSol
                print '{:<4}    {:7.2f}    {:7.2f}    {:4.2f}'.format(str(vesselId), volumeInit, volumeSol, diff)
            
            print ' ----------------------------------- '
            print 'Boundarys        in (ml)   out (ml)  '
            totalIn = 0
            totalOut = 0            
            for boundaryList in self.boundarys.itervalues():    
                for boundary in boundaryList:
                    #print '  volume in   (ml)  {:.2f}'.format((abs(boundary.BloodVolumen[0])*1.e6))
                    #print '  volume out  (ml)  {:.2f}'.format((abs(boundary.BloodVolumen[1])*1.e6))
                    #print '  volume diff (ml)  {:.2f}'.format(abs(abs(boundary.BloodVolumen[0])*1.e6 - abs(boundary.BloodVolumen[1])*1.e6))
                    print '{:<12}     {:6.2f}    {:4.2f}'.format(boundary.name,
                                                                 abs(boundary.BloodVolumen[0])*1.e6, 
                                                                 abs(boundary.BloodVolumen[1])*1.e6 ) 
                                                             #abs(abs(boundary.BloodVolumen[0])*1.e6 - abs(boundary.BloodVolumen[1])*1.e6))
                    #print '   {} '.format(boundary.type)
                    totalIn  += abs(boundary.BloodVolumen[0])*1.e6
                    totalOut += abs(boundary.BloodVolumen[1])*1.e6
            print ' ----------------------------------- '
            print 'total'
            print '  vessels initial  (ml)  {:4.2f}'.format(vesselsInit)  
            print '  vessels solution (ml)  {:4.2f}'.format(vesselsSol)  
            print '  boundarys in     (ml)  {:4.2f}'.format(totalIn)
            print '  boundarys out    (ml)  {:4.2f}'.format(totalOut)
            print '                    ------------'
            print '  volume diff      (ml)  {:.2f}'.format(vesselsInit - vesselsSol + totalIn - totalOut)
        
        # postprocessing: calculate wave speed
        for vessel in self.vessels.itervalues():
            vessel.postProcessing()  
            
        ## stop realtime visualisation
        for communicator in self.communicators.itervalues():           
            try: communicator.stopRealtimeViz()
            except: pass
            
        ### garbage collection
        gc.collect()
        
        del self.numericalObjects
        del self.fields
        del self.connections
        del self.boundarys
        del self.communicators
        
        
        ## create Mimic of the old solutionDataSets from solution data saved in each vessel instance
        # this will be removed soon!!!
        
         # solution of the system over time {vesselID: [ [solution at N nodes]<-one array for each timePoint , ...  ]  }
        # old and will be removed soon
        P  = {}
        Q  = {}
        A  = {}
        c  = {}
        for vesselId,vessel in self.vessels.iteritems():
            P[vesselId]    = self.vessels[vesselId].Psol
            Q[vesselId]    = self.vessels[vesselId].Qsol
            A[vesselId]    = self.vessels[vesselId].Asol
            c[vesselId]    = self.vessels[vesselId].csol
        
        return P,Q,A,c