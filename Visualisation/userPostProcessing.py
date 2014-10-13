import sys,os
cur = os.path.dirname( os.path.realpath( __file__ ) )

sys.path.append(cur+'/../UtilityLib')
from moduleXML import loadNetworkFromXML 
from moduleStartUp import parseOptions
from modulePickle import loadSolutionDataFile


sys.path.append(cur+'/../NetworkLib')
from classVascularNetwork import VascularNetwork 

from optparse import OptionParser
import cPickle
import numpy as np

from pprint import pprint as pp

__author__ = "Vinzenz Gregor Eck"
__version__ = "0.3"


### file to help user post processing

### 1. open solution data file after a simulation

## use -f and -n to define solution file you want to open, or use the inbulid menu

optionsDict = parseOptions(['f','n','c'], visualisationOnly = True)
    
networkName           = optionsDict['networkName']
dataSetNumber         = optionsDict['dataNumber']

##  open data file choosed above
try:
    print " Try to open network {} with data number {}".format(networkName, dataSetNumber)
    vascularNetwork, solutionDataSets, simulationCaseDescriptions = loadSolutionDataFile(networkName, [dataSetNumber])
    solutionData = solutionDataSets[0]
except:
    print "Error could not open solution data with data number {} of network {}".format(dataSetNumber,networkName)
    exit()
    

### 2. available data
print "\n 2. available data \n"
## 2.1 vascularNetwork  
## vascularNetwork is an instance of the vascularNetwork class containting all data which was available for the simulation
## as e.g. totalTime, dt, all vessel instances etc.
## see also in classVascularNetwork for methods and variables available

## Usage:

print "timeStep of simulation: vascularNetwork.dt", vascularNetwork.dt
try:
    print "gridpoints of vessel with Id 0: vascularNetwork.vessels[0].N", vascularNetwork.vessels[0].N
except: pass

## 2.1 solutionData
## is a nested dictionary containting the solution of the simulation
## the first rank contains['Pressure', 'Area', 'Flow', 'Name', 'WaveSpeed']
## see as well:
print solutionData.keys()
## the second rank contains a dict for all vesselIds 
## e.g.
vesselId = 0
solutionType = "Pressure"
print solutionData[solutionType][vesselId]
## in the third rank, a numpy matrix with solutions for each timestep and each grid node are saved:
## solution for all grid nodes at a point in time:
print solutionData[solutionType][vesselId][0]
## solution for all points in time at grid node 0:
print solutionData[solutionType][vesselId][:,[0]]


