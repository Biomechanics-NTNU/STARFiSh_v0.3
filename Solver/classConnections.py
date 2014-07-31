import numpy as np
from numpy.linalg import solve
from scipy.optimize import fsolve

from pprint import pprint as pp       

from copy import copy as copy       

__author__ = "Vinzenz Gregor Eck"
__version__ = "0.3"

class Link():        
    '''
    Link object represends the connection between two vessels
    '''
    def __init__(self, mother, motherSys, 
                     daughter, daughterSys,
                     n, dt, rigidAreas, solvingScheme):
        self.type = 'Link'
        
        self.name = ' '.join(['Link',str(mother.Id),str(daughter.Id)])
        
        self.rho             = []
        self.systemEquations = []
        self.z               = []
        self.A_func          = []
        self.positions       = []
        self.names           = []
#       self.vz = []

        self.dt = dt
        self.n  = n
        
        # equations to solve in f solve
        self.fsolveFunction = None
        self.jacobiMatrix = None
        
        #initialize Vessels
        # mother vessel
        self.rho.append(mother.rho)
        self.z.append(mother.z)
        self.systemEquations.append(motherSys)
        self.positions.append(-1)
        self.names.append(mother.Id)
        self.A_func.append(mother.A_nID)
#        self.vz.append(-1)
        # SolutionVariables
        self.P_mother = mother.Psol
        self.Q_mother = mother.Qsol
        self.A_mother = mother.Asol
               
        # daughter vessel
        self.rho.append(daughter.rho)
        self.z.append(daughter.z)
        self.systemEquations.append(daughterSys)
        self.positions.append(0)
        self.names.append(daughter.Id)
        self.A_func.append(daughter.A_nID)
#        self.vz.append(1)
        # SolutionVariables
        self.P_daughter = daughter.Psol
        self.Q_daughter = daughter.Qsol
        self.A_daughter = daughter.Asol
            
#        if   rigidAreas == '02': 
#            self.fsolveFunction = self.fsolveConnectionSys0
##            self.jacobiMatrix = self.jacobiMatrixSys0
#        elif rigidAreas == '2':  
#            self.fsolveFunction = self.fsolveConnectionSys1
#            self.jacobiMatrix = self.jacobiMatrixSys1
#       else: 
#            print"ERROR: classConnections: EquSys not properly defined! system exit"
#            exit()

        self.rigidAreas = rigidAreas
        
        # Define the call function depending on the solving Scheme
        if solvingScheme == "Linear": 
            self.__call__ = self.callLinear
        else:
            print "ERROR Connections wrong solving scheme!", solvingScheme; exit()
    
        ## benchamark Test variables
        self.nonLin = False
        self.updateL = False
        self.sumQErrorCount = 0
        self.maxQError = 0
        self.maxPErrorNonLin = 0 
        self.maxPError = 0
        self.sumPErrorCount = 0
        self.sumPErrorNonLinCount = 0
    
    def callLinear(self):
        '''
        Call function for vessel-vessel connection
        '''        
        dt = self.dt
        n = self.n[0]
        pos1 = self.positions[0]
        pos2 = self.positions[1]
        
        P1 = self.P_mother[n]
        Q1 = self.Q_mother[n]
        A1 = self.A_mother[n]
        
        P2 = self.P_daughter[n]
        Q2 = self.Q_daughter[n]
        A2 = self.A_daughter[n]
        
        P1o = P1[pos1]
        Q1o = Q1[pos1]
        P2o = P2[pos2]
        Q2o = Q2[pos2]
        
        # update system equation and store L1
        self.systemEquations[0].updateLARL(P1,Q1,A1,idArray=[pos1]) 
        #L1 = self.systemEquations[0].L[pos1][pos1+1]
        R1 = self.systemEquations[0].R[pos1]
        domega1_1 = self.systemEquations[0].domega[pos1]
        
        self.systemEquations[1].updateLARL(P2,Q2,A2,idArray=[pos2]) 
        #L2 = self.systemEquations[1].L[pos2][pos2+1]
        R2 = self.systemEquations[1].R[pos2]
        domega2_2 = self.systemEquations[1].domega[pos2]
    
            
        # local R matrices
        R1_11 = R1[0][0]
        R1_12 = R1[0][1]
        R1_21 = R1[1][0]
        R1_22 = R1[1][1]
        
        R2_11 = R2[0][0]
        R2_12 = R2[0][1]
        R2_21 = R2[1][0]
        R2_22 = R2[1][1]
                
        ###### Linear approach
        denom = (R1_12*R2_21-R1_22*R2_11)
        # -1 reflectionCoeff mother->daugther
        alpha1 = -(R1_11*R2_21-R1_21*R2_11)/denom
        # +1 transmission daugther->mother
        alpha2 = -(R2_11*R2_22-R2_21*R2_12)/denom 
        # +1 transmission daugther->mother
        beta1 = -(R1_11*R1_22-R1_12*R1_21)/denom 
        # -1 reflectionCoeff daugther->mother
        beta2 = -(R1_12*R2_22-R2_12*R1_22)/denom
        
        #print 'cC153 alphas',alpha1,alpha2,beta1,beta2
                
        domega1_2 = alpha1 * domega1_1 + alpha2 * domega2_2
        domega2_1 = beta1  * domega1_1 + beta2  * domega2_2
                
        P1_new = P1o + R1_11*domega1_1 + R1_12*domega1_2
        Q1_new = Q1o + R1_21*domega1_1 + R1_22*domega1_2
    
        P2_new = P2o + R2_11*domega2_1 + R2_12*domega2_2
        Q2_new = Q2o + R2_21*domega2_1 + R2_22*domega2_2
        
        # apply calculated values to next time step
        self.P_mother[n+1][pos1]   = P1_new
        self.Q_mother[n+1][pos1]   = Q1_new
        self.P_daughter[n+1][pos2] = P2_new
        self.Q_daughter[n+1][pos2] = Q2_new
        
        if P1_new < 0 or P2_new < 0:
            print "ERROR: Connection: {} calculated negativ pressure at time {} (n {},dt {}), exit system".format(self.name,n*dt,n,dt)
            print P1_new, P2_new
            exit()
                
        # calculate new areas
        if self.rigidAreas == False:
            A1n = self.A_func[0]([P1_new],pos1)
            A2n = self.A_func[1]([P2_new],pos2)          
        else:
            A1n = A1[pos1]
            A2n = A2[pos2]
        # apply Areas
        self.A_mother[n+1][pos1]   = A1n
        self.A_daughter[n+1][pos2] = A2n
              
        # Error estimation
        try: sumQError = abs(abs(Q1_new)-abs(Q2_new))/abs(Q1_new)
        except: sumQError = 0.0
        if sumQError > 0.0: 
            self.sumQErrorCount = self.sumQErrorCount+1
        if sumQError > self.maxQError:
            self.maxQError  = sumQError
        #print self.name, ' Error cons mass',  sumQError, self.maxQError ,' - ', n, self.sumQErrorCount
        
        # Linear presssure equation
        sumPError = abs(abs(P1_new)-abs(P2_new))/abs(P1_new)
        if sumPError > 0.0: 
            self.sumPErrorCount = self.sumPErrorCount+1
        if sumPError > self.maxPError:
            self.maxPError  = sumPError
        if sumQError > 1.e-5:
            print "ERROR: Connection: {} to high error in conservation of mass at time {} (n {},dt {}), exit system".format(self.name,n*dt,n,dt)
            print sumQError
            exit()
                                               
        ## Non linear pressure equation
        sumPErrorNonLin = abs(abs(P1_new+1000*0.5*(Q1_new/A1n)**2)-abs(abs(P2_new+1000*0.5*(Q2_new/A2n)**2)))/abs(P1_new+0.5*(Q1_new/A1n)**2)
        if sumPErrorNonLin > 0.0: 
            self.sumPErrorNonLinCount = self.sumPErrorNonLinCount+1
        if sumPErrorNonLin > self.maxPErrorNonLin:
            self.maxPErrorNonLin  = sumPErrorNonLin
        if sumPError > 1.e-10:
            print "ERROR: Connection: {} to high error in conservation of pressure at time {} (n {},dt {}), exit system".format(self.name,n*dt,n,dt)
            print sumPError
            exit()    
    
class Bifurcation():
    
    def __init__(self, mother, motherSys,
                       leftDaughter, leftDaughterSys,
                       rightDaughter, rightDaughterSys, 
                       n, dt, rigidAreas, solvingScheme):
        # vessel variables initially set, constant through simulation
        self.type = 'Bifurcation'
        
        self.name = ' '.join(['Bifurcation',str(mother.Id),str(leftDaughter.Id),str(rightDaughter.Id)])
        
        #System Variables
        self.dt = dt
        self.n = n
        
        self.rho = []
        self.systemEquations = []
        self.z = []
        self.A_func = []
        self.positions =[]
#         self.vz = []
        self.names = []
        
#         # equations to solve in f solve
#         self.fsolveFunction = None
#         self.jacobiMatrix = None
        
        ###initialize
        ##mother branch
        self.rho.append(mother.rho)
        self.z.append(mother.z)
        self.systemEquations.append(motherSys)
        self.positions.append(-1)
        self.names.append(mother.Id)
        self.A_func.append(mother.A_nID)
#        self.vz.append(-1)
        #SolutionVariables
        self.P_leftMother = mother.Psol
        self.Q_leftMother = mother.Qsol
        self.A_leftMother = mother.Asol
              
        
        ##left daughter
        self.rho.append(leftDaughter.rho)
        self.z.append(leftDaughter.z)
        self.systemEquations.append(leftDaughterSys)
        self.positions.append(0)
        self.names.append(leftDaughter.Id)
        self.A_func.append(leftDaughter.A_nID)
#        self.vz.append(1)
        #SolutionVariables
        self.P_leftDaughter = leftDaughter.Psol
        self.Q_leftDaughter = leftDaughter.Qsol
        self.A_leftDaughter = leftDaughter.Asol
        
        ##right daughter
        self.rho.append(rightDaughter.rho)
        self.z.append(rightDaughter.z)
        self.systemEquations.append(rightDaughterSys)
        self.positions.append(0)
        self.names.append(rightDaughter.Id)
        self.A_func.append(rightDaughter.A_nID)
#        self.vz.append(1)
        #SolutionVariables
        self.P_rightDaughter = rightDaughter.Psol
        self.Q_rightDaughter = rightDaughter.Qsol
        self.A_rightDaughter = rightDaughter.Asol
        
        self.rigidAreas = rigidAreas
    
        # Define the call function depending on the solving Scheme
        if solvingScheme == "Linear": 
            self.__call__ = self.callLinear
        else:
            print "ERROR Connections wrong solving scheme!", solvingScheme; exit()
        
        ## benchamark Test variables
        self.sumQErrorCount = 0
        self.maxQError = 0
        self.maxPErrorNonLin = 0 
        self.maxPError = 0
        self.sumPErrorCount = 0
        self.sumPErrorNonLinCount = 0
    
    def callLinear(self):
        '''
        Call function for vessel-vessel connection
        '''        
        dt = self.dt
        n = self.n[0]
        pos1 = self.positions[0]
        pos2 = self.positions[1]
        pos3 = self.positions[2]
        
        P1 = self.P_leftMother[n]
        Q1 = self.Q_leftMother[n]
        A1 = self.A_leftMother[n]
        
        P2 = self.P_leftDaughter[n]
        Q2 = self.Q_leftDaughter[n]
        A2 = self.A_leftDaughter[n]
        
        P3 = self.P_rightDaughter[n]
        Q3 = self.Q_rightDaughter[n]
        A3 = self.A_rightDaughter[n]
        
        P1o = P1[pos1]
        Q1o = Q1[pos1]
        P2o = P2[pos2]
        Q2o = Q2[pos2]
        P3o = P3[pos3]
        Q3o = Q3[pos3]
        
        ## update LARL
        self.systemEquations[0].updateLARL(P1,Q1,A1,idArray=[pos1])
        L1 = self.systemEquations[0].L[pos1][pos1+1]
        R1 = self.systemEquations[0].R[pos1]
        domega1_1 = self.systemEquations[0].domega[pos1]
        
        self.systemEquations[1].updateLARL(P2,Q2,A2,idArray=[pos2])
        L2 = self.systemEquations[1].L[pos2][pos2+1]
        R2 = self.systemEquations[1].R[pos2]
        domega2_2 = self.systemEquations[1].domega[pos2]
        
        self.systemEquations[2].updateLARL(P3,Q3,A3,idArray=[pos3])                                    
        L3 = self.systemEquations[2].L[pos3][pos3+1]
        R3 = self.systemEquations[2].R[pos2]
        domega3_2 = self.systemEquations[2].domega[pos3]
        
        # local R matrices
        R1_11 = R1[0][0]
        R1_12 = R1[0][1]
        R1_21 = R1[1][0]
        R1_22 = R1[1][1]
        
        R2_11 = R2[0][0]
        R2_12 = R2[0][1]
        R2_21 = R2[1][0]
        R2_22 = R2[1][1]
        
        R3_11 = R3[0][0]
        R3_12 = R3[0][1]
        R3_21 = R3[1][0]
        R3_22 = R3[1][1]
        
        ###### Linear approach
        denom = R1_12 * R2_11 * R3_21 + R1_12 * R2_21 * R3_11 - R1_22 * R2_11 * R3_11 
        
        alpha1 = -(R1_11 * R2_11 * R3_21 + R1_11 * R2_21 * R3_11 - R1_21 * R2_11 * R3_11)/denom
        alpha2 = -(R2_11 * R2_22 * R3_11 - R2_12 * R2_21 * R3_11)/denom
        alpha3 = -(R2_11 * R3_22 * R3_11 - R2_11 * R3_12 * R3_21)/denom
        
        beta1 = -(R1_11 * R1_22 * R3_11 - R1_12 * R1_21 * R3_11)/denom
        beta2 = -(R1_12 * R2_12 * R3_21 + R1_12 * R2_22 * R3_11 - R1_22 * R2_12 * R3_11)/denom
        beta3 = -(R1_12 * R3_11 * R3_22 - R1_12 * R3_12 * R3_21)/denom
        
        gamma1 = -(R1_11 * R1_22 * R2_11 - R1_12 * R1_21 * R2_11)/denom
        gamma2 = -(R1_12 * R2_11 * R2_22 - R1_12 * R2_12 * R2_21)/denom
        gamma3 = -(R1_12 * R2_11 * R3_22 + R1_12 * R2_21 * R3_12 - R1_22 * R2_11 * R3_12)/denom
        
        domega1_2 = alpha1 * domega1_1  + alpha2 * domega2_2 + alpha3 * domega3_2 
        domega2_1 = beta1  * domega1_1  + beta2  * domega2_2 + beta3  * domega3_2 
        domega3_1 = gamma1 * domega1_1  + gamma2 * domega2_2 + gamma3 * domega3_2 
        
        P1_new = P1o + (R1_11*domega1_1 + R1_12*domega1_2)
        Q1_new = Q1o + (R1_21*domega1_1 + R1_22*domega1_2)
    
        P2_new = P2o + (R2_11*domega2_1 + R2_12*domega2_2)
        Q2_new = Q2o + (R2_21*domega2_1 + R2_22*domega2_2)
        
        P3_new = P3o + (R3_11*domega3_1 + R3_12*domega3_2)
        Q3_new = Q3o + (R3_21*domega3_1 + R3_22*domega3_2)
               
        
        # apply calculated values to next time step
        self.P_leftMother[n+1][pos1]    = P1_new
        self.Q_leftMother[n+1][pos1]    = Q1_new
        self.P_leftDaughter[n+1][pos2]  = P2_new
        self.Q_leftDaughter[n+1][pos2]  = Q2_new
        self.P_rightDaughter[n+1][pos3] = P3_new
        self.Q_rightDaughter[n+1][pos3] = Q3_new
        
        if P1_new < 0 or P2_new < 0 or P3_new < 0:
            print "ERROR: Connection: {} calculated negativ pressure at time {} (n {},dt {}), exit system".format(self.name,n*dt,n,dt)
            print P1_new, P2_new, P3_new
            exit()
        
        
        # calculate new areas
        if self.rigidAreas == False:
            A1n = self.A_func[0]([P1_new],pos1)
            A2n = self.A_func[1]([P2_new],pos2)
            A3n = self.A_func[2]([P3_new],pos3)
        else:
            A1n = A1[pos1]
            A2n = A2[pos2]
            A3n = A3[pos3] 
           
        self.A_leftMother[n+1][pos1]       = A1n
        self.A_leftDaughter[n+1][pos2]     = A2n       
        self.A_rightDaughter[n+1][pos3]    = A3n  
            
        ## non linear error        
        try: sumQError = abs(Q1_new-Q2_new-Q3_new)/abs(Q1_new)
        except: sumQError = 0.0
        if sumQError > 0.0: 
            self.sumQErrorCount = self.sumQErrorCount+1
        if sumQError > self.maxQError:
            self.maxQError  = sumQError
        #print self.name,' \n Error cons mass',  sumQError, self.maxQError ,' - ', n, self.sumQErrorCount
        if sumQError > 1.e-5:
            print "ERROR: Connection: {} to high error in conservation of mass at time {} (n {},dt {}), exit system".format(self.name,n*dt,n,dt)
            print sumQError
            exit()
        
        sumPError = abs(P1_new-P2_new)/abs(P1_new)
        if sumPError > 0.0: 
            self.sumPErrorCount = self.sumPErrorCount+1
        if sumPError > self.maxPError:
            self.maxPError  = sumPError
        #print self.name,' Error P lin    ',  sumPError, self.maxPError ,' - ', n, self.sumPErrorCount
        if sumPError > 1.e-10:
            print "ERROR: Connection: {} to high error in conservation of pressure at time {} (n {},dt {}), exit system".format(self.name,n*dt,n,dt)
            print sumPError
            exit()
        
        sumPErrorNonLin = abs(P1_new+500*(Q1_new/A1n)**2-(P2_new+500*(Q2_new/A2n)**2))/abs(P1_new+0.5*(Q1_new/A1n)**2)
        if sumPErrorNonLin > 0.0: 
            self.sumPErrorNonLinCount = self.sumPErrorNonLinCount+1
        if sumPErrorNonLin > self.maxPErrorNonLin:
            self.maxPErrorNonLin  = sumPErrorNonLin
        #print self.name,' Error P non lin',  sumPErrorNonLin, self.maxPErrorNonLin ,' - ', n, self.sumPErrorNonLinCount
        
                
class Anastomosis():
    
    def __init__(self, leftMother, leftMotherSys,
                       rightMother, rightMotherSys,
                       daughter, daughterSys, 
                       n, dt, rigidAreas, solvingScheme):
        
        
        # vessel variables initially set, constant through simulation
        self.type = 'Anastomosis'
        
        self.name = ' '.join(['Anastomosis',str(leftMother.Id),str(rightMother.Id),str(daughter.Id)])
        
        #System Variables
        self.dt = dt
        self.n = n
        
        self.rho = []
        self.systemEquations = []
        self.z = []
        self.A_func = []
        self.positions =[]
#         self.vz = []
        self.names = []
        
#         # equations to solve in f solve
#         self.fsolveFunction = None
#         self.jacobiMatrix = None
        
        ###initialize
        ##left mother branch
        self.rho.append(leftMother.rho)
        self.z.append(leftMother.z)
        self.systemEquations.append(leftMotherSys)
        self.positions.append(-1)
        self.names.append(leftMother.Id)
        self.A_func.append(leftMother.A_nID)
        #SolutionVariables
        self.P_leftMother = leftMother.Psol
        self.Q_leftMother = leftMother.Qsol
        self.A_leftMother = leftMother.Asol
        
        ##left mother branch
        self.rho.append(rightMother.rho)
        self.z.append(rightMother.z)
        self.systemEquations.append(rightMotherSys)
        self.positions.append(-1)
        self.names.append(rightMother.Id)
        self.A_func.append(rightMother.A_nID)
        #SolutionVariables
        self.P_rightMother = rightMother.Psol
        self.Q_rightMother = rightMother.Qsol
        self.A_rightMother = rightMother.Asol
        
        ##right daughter branch
        self.rho.append(daughter.rho)
        self.z.append(daughter.z)
        self.systemEquations.append(daughterSys)
        self.positions.append(0)
        self.names.append(daughter.Id)
        self.A_func.append(daughter.A_nID)
#        self.vz.append(1)
        #SolutionVariables
        self.P_daughter = daughter.Psol
        self.Q_daughter = daughter.Qsol
        self.A_daughter = daughter.Asol

        self.rigidAreas = rigidAreas
    
        # Define the call function depending on the solving Scheme
        if solvingScheme == "Linear": 
            self.__call__ = self.callLinear
        else:
            print "ERROR Connections wrong solving scheme!", solvingScheme; exit()
        
        ## benchamark Test variables
        self.sumQErrorCount = 0
        self.maxQError = 0
        self.maxPErrorNonLin = 0 
        self.maxPError = 0
        self.sumPErrorCount = 0
        self.sumPErrorNonLinCount = 0
    
    def callLinear(self):
        '''
        Call function for vessel-vessel connection
        '''        
        dt = self.dt
        n = self.n[0]
        pos1 = self.positions[0]
        pos2 = self.positions[1]
        pos3 = self.positions[2]
        
        P1 = self.P_leftMother[n]
        Q1 = self.Q_leftMother[n]
        A1 = self.A_leftMother[n]
        
        P2 = self.P_rightMother[n]
        Q2 = self.Q_rightMother[n]
        A2 = self.A_rightMother[n]
        
        P3 = self.P_daughter[n]
        Q3 = self.Q_daughter[n]
        A3 = self.A_daughter[n]
        
        P1o = P1[pos1]
        Q1o = Q1[pos1]
        P2o = P2[pos2]
        Q2o = Q2[pos2]
        P3o = P3[pos3]
        Q3o = Q3[pos3]
                           
                
        # update LARL
        self.systemEquations[0].updateLARL(P1,Q1,A1,idArray=[pos1]) #
        L1 = self.systemEquations[0].L[pos1][pos1+1]
        R1 = self.systemEquations[0].R[pos1]
        domega1_1 = self.systemEquations[0].domega[pos1]
        
        self.systemEquations[1].updateLARL(P2,Q2,A2,idArray=[pos2])
        L2 = self.systemEquations[1].L[pos2][pos2+1]
        R2 = self.systemEquations[1].R[pos2]
        domega2_1 = self.systemEquations[1].domega[pos2]
        
        self.systemEquations[2].updateLARL(P3,Q3,A3,idArray=[pos3])
        L3 = self.systemEquations[2].L[pos3][pos3+1]
        R3 = self.systemEquations[2].R[pos2]
        domega3_2 = self.systemEquations[2].domega[pos3]
            
        # local R matrices
        R1_11 = R1[0][0]
        R1_12 = R1[0][1]
        R1_21 = R1[1][0]
        R1_22 = R1[1][1]
        
        R2_11 = R2[0][0]
        R2_12 = R2[0][1]
        R2_21 = R2[1][0]
        R2_22 = R2[1][1]
        
        R3_11 = R3[0][0]
        R3_12 = R3[0][1]
        R3_21 = R3[1][0]
        R3_22 = R3[1][1]
        
        ###### Linear approach
        
        ####### change?!
        denom = R1_12 * R2_12 * R3_21 - R1_12 * R2_22 * R3_11 - R1_22 * R2_12 * R3_11 
        
        alpha1 = -( R1_11 * R2_12 * R3_21 - R1_11 * R2_22 * R3_11 - R1_21 * R2_12 * R3_11)/denom
        alpha2 = -( R2_11 * R2_22 * R3_11 - R2_12 * R2_21 * R3_11)/denom
        alpha3 = -( R2_12 * R3_22 * R3_11 - R2_12 * R3_12 * R3_21)/denom
        
        beta1 = -(R1_11 * R1_22 * R3_11 - R1_12 * R1_21 * R3_11)/denom
        beta2 = -(R1_12 * R2_11 * R3_21 - R1_12 * R2_21 * R3_11 - R1_22 * R2_11 * R3_11)/denom
        beta3 = -(R1_12 * R3_11 * R3_22 - R1_12 * R3_12 * R3_21)/denom
        
        gamma1 = -( R1_11 * R1_22 * R2_12 - R1_12 * R1_21 * R2_12)/denom
        gamma2 = -( R1_12 * R2_11 * R2_22 + R1_12 * R2_12 * R2_21)/denom
        gamma3 = -( R1_12 * R2_12 * R3_22 - R1_12 * R2_22 * R3_12 - R1_22 * R2_12 * R3_12)/denom
        
        ############        
        domega1_2 = (alpha1 * domega1_1  + alpha2 * domega2_1 + alpha3 * domega3_2 )
        domega2_2 = (beta1  * domega1_1  + beta2  * domega2_1 + beta3  * domega3_2 )
        domega3_1 = (gamma1 * domega1_1  + gamma2 * domega2_1 + gamma3 * domega3_2 )
        
        P1_new = P1o + (R1_11*domega1_1 + R1_12*domega1_2)
        Q1_new = Q1o + (R1_21*domega1_1 + R1_22*domega1_2)
    
        P2_new = P2o + (R2_11*domega2_1 + R2_12*domega2_2)
        Q2_new = Q2o + (R2_21*domega2_1 + R2_22*domega2_2)
        
        P3_new = P3o + (R3_11*domega3_1 + R3_12*domega3_2)
        Q3_new = Q3o + (R3_21*domega3_1 + R3_22*domega3_2)
               
        
        # apply calculated values to next time step solution
        self.P_leftMother[n+1][pos1]  = P1_new
        self.Q_leftMother[n+1][pos1]  = Q1_new
        self.P_rightMother[n+1][pos2] = P2_new
        self.Q_rightMother[n+1][pos2] = Q2_new
        self.P_daughter[n+1][pos3]    = P3_new
        self.Q_daughter[n+1][pos3]    = Q3_new
        
        if P1_new < 0 or P2_new < 0 or P3_new < 0:
            print "ERROR: Connection: {} calculated negativ pressure at time {} (n {},dt {}), exit system".format(self.name,n*dt,n,dt)
            print P1_new, P2_new, P3_new
            exit()
                
        # calculate new areas
        if self.rigidAreas == False:
            A1n = self.A_func[0]([P1_new],pos1)
            A2n = self.A_func[1]([P2_new],pos2)
            A3n = self.A_func[2]([P3_new],pos3)
        else:
            A1n = A1[pos1]
            A2n = A2[pos2]
            A3n = A3[pos3] 
           
        self.A_leftMother[n+1][pos1]   = A1n
        self.A_rightMother[n+1][pos2]  = A2n       
        self.A_daughter[n+1][pos3]     = A3n  
            
        ## non linear error        
        try: sumQError = abs(Q1_new+Q2_new-Q3_new)#/abs(Q1_new)
        except: sumQError = 0.0
        if sumQError > 0.0: 
            self.sumQErrorCount = self.sumQErrorCount+1
        if sumQError > self.maxQError:
            self.maxQError  = sumQError
        #print self.name,' \n Error cons mass',  sumQError, self.maxQError ,' - ', n, self.sumQErrorCount
        if sumQError > 1.e-5:
            print "ERROR: Connection: {} to high error in conservation of mass at time {} (n {},dt {}), exit system".format(self.name,n*dt,n,dt)
            print sumQError, ' <- Q1,Q2,Q3 ',Q1_new, Q2_new, Q3_new
            #exit()
        
        sumPError = abs(P1_new-P3_new)/abs(P3_new)
        if sumPError > 0.0: 
            self.sumPErrorCount = self.sumPErrorCount+1
        if sumPError > self.maxPError:
            self.maxPError  = sumPError
        #print self.name,' Error P lin    ',  sumPError, self.maxPError ,' - ', n, self.sumPErrorCount
        if sumPError > 1.e-5:
            print "ERROR: Connection: {} to high error in conservation of pressure at time {} (n {},dt {}), exit system".format(self.name,n*dt,n,dt)
            print sumPError, ' <- P1,P2,P3 ',P1_new, P2_new, P3_new
            #exit()
        
        sumPErrorNonLin = 1050./2.*(Q1_new/A1n)**2#abs(P1_new+500*(Q1_new/A1n)**2-(P2_new+500*(Q2_new/A2n)**2))/abs(P1_new+0.5*(Q1_new/A1n)**2)
        if sumPErrorNonLin > 0.0: 
            self.sumPErrorNonLinCount = self.sumPErrorNonLinCount+1
        if sumPErrorNonLin > self.maxPErrorNonLin:
            self.maxPErrorNonLin  = sumPErrorNonLin
        print self.name,' Error P non lin',  sumPErrorNonLin, self.maxPErrorNonLin ,' - ', n, self.sumPErrorNonLinCount
        
        print self.name,'dynamic Pressures',1050./2.*(Q1_new/A1n)**2,1050./2.*(Q2_new/A2n)**2,1050./2.*(Q3_new/A3n)**2, '--',1050./2.*(Q1_new/A1n)**2+-1050./2.*(Q3_new/A3n)**2
        
        