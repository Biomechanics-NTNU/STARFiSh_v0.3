########################################################################################
#                            STARFiSh v0.3
########################################################################################
## 
# created by Vinzenz Eck vinzenz.eck@mytum.de
##
#---------------------------------------------------------------------------------------#
import time 
import sys,os
# set the path relative to THIS file not the executing file!
cur = os.path.dirname( os.path.realpath( __file__ ) )
#sys.path.append(cur+'/Solver')

sys.path.append(cur+'/NetworkLib')
from classVascularNetwork import VascularNetwork

sys.path.append(cur+'/Solver')
from class1DflowSolver import FlowSolver

sys.path.append(cur+'/UtilityLib')
from moduleXML import loadNetworkFromXML
from moduleStartUp import parseOptions
from modulePickle import saveSolutionDataFile
from modulePickle import updateSimulationDescriptions
from modulePickle import loadSolutionDataFile

sys.path.append(cur+'/Visualisation')
from class3dVisualisation import Visualisation3D

import pprint
import matplotlib.pyplot as plt

import cPickle

import gc

import subprocess


__author__ = "Vinzenz Gregor Eck"
__version__ = "0.3"

def main():
    print ""
    print '====================================='
    print '#          STARFiSh_v0.3            #'
    print '====================================='
    
    optionsDict = parseOptions(['f','n','d','s','v','r'])
    
    networkName           = optionsDict['networkName']
    save                  = optionsDict['save']
    dataNumber            = optionsDict['dataNumber']
    simulationDescription = optionsDict['simulationDescription']
    vizOutput             = optionsDict['vizOutput']
    resimulate            = optionsDict['resimulate']
    
    filename = str(networkName+'.xml')
        
    print '____________Simulation_______________'
    print '%-20s %s' % ('Network name',networkName)
    print '%-20s %s' % ('Data number', dataNumber)
    print '%-20s %s' % ('Save simulation', save)
    print '%-20s %s' % ('Case description', simulationDescription)
    print '%-20s %s' % ('Resimulate', resimulate)
    print '%-20s %s' % ('Visualisationmode', vizOutput)
    
    
    # load network from the path!
    if resimulate == False:
        vascularNetwork = loadNetworkFromXML(filename=filename)
    else:
        vascularNetwork,solutionDataSets,simulationCaseDescriptions = loadSolutionDataFile(networkName, [dataNumber])
        simulationDescription = simulationCaseDescriptions[0][0]
        print simulationCaseDescriptions
        print simulationDescription,type(simulationDescription)
    
    if vascularNetwork == None: exit()
    vascularNetwork.evaluateConnections()
    
    timeSolverInitStart = time.clock()
    #initialize Solver
    flowSolver = FlowSolver(vascularNetwork)
    timeSolverInit = time.clock()-timeSolverInitStart
    timeSolverSolveStart = time.clock()
    #solve the system
    P,Q,A,c = flowSolver.solve()
    timeSolverSolve = time.clock()-timeSolverSolveStart
    
    minutesInit = int(timeSolverInit/60.)
    secsInit = timeSolverInit-minutesInit*60.
    minutesSolve = int(timeSolverSolve/60.)
    secsSolve = timeSolverSolve-minutesSolve*60.
    
    
    #print '====================================='
    print '____________ Solver time _____________'
    print 'Initialisation: {} min {} sec'.format(minutesInit,secsInit)
    print 'Solving:        {} min {} sec'.format(minutesSolve,secsSolve)
    print '====================================='
    
          
    solutionDataSave = [{'Pressure': P, 'Flow': Q, 'Area': A, 'WaveSpeed':c, 'Name': str('simulation_'+dataNumber) }]
    networkName = filename.split('.')[0]

    saveSolutionDataFile(networkName,dataNumber,vascularNetwork,solutionDataSave,simulationDescription)
    
    del flowSolver
    gc.collect()
    
    updateSimulationDescriptions(networkName,dataNumber)
    
    gc.collect()
    
    if vizOutput == "2D":
        string = ' '.join(['python',cur+'/Visualisation/class2dVisualisation.py','-f',vascularNetwork.name, '-n',str(dataNumber)])                
        subprocess.Popen(string, shell=True)
        
    
    if vizOutput == "3D":
        string = ' '.join(['python',cur+'/Visualisation/class3dVisualisation.py','-f',vascularNetwork.name, '-n',str(dataNumber), '-c True']) 
        subprocess.Popen(string, shell=True)
        
        
    if vizOutput == "2D+3D":
           
        string1 = ' '.join(['python',cur+'/Visualisation/class2dVisualisation.py','-f',vascularNetwork.name, '-n',str(dataNumber), '-c True']) 
        string2 = ' '.join(['python',cur+'/Visualisation/class3dVisualisation.py','-f',vascularNetwork.name, '-n',str(dataNumber), '-c True']) 
        
        viz2d = subprocess.Popen(string1, shell = True )
        viz3d = subprocess.Popen(string2, shell = True )
        
        while True:
            
            if viz2d.poll() != None:
                viz3d.terminate()
                exit()
                
            if viz3d.poll() != None:
                viz2d.terminate()
                exit()
        
if __name__ == '__main__':
    
    profiling = False
    callGraph = False
    
    if profiling == True:
        import cProfile
        cProfile.run('main()')
        
    elif callGraph == True:
        import pycallgraph
        pycallgraph.start_trace()
        main()
        pycallgraph.make_dot_graph('vascular1DFlowCallGraph.png')
    else:
        main()
