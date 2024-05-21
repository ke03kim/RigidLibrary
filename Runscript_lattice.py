import random
import sys, os, glob

# Load pebble game and dynamical matrix libraries
import Configuration as CF
import Pebbles as PB
import Hessian as HS
import Analysis as AN
import Tiling as TY

import matplotlib.pyplot as plt
import numpy as np
import pickle
import pandas as pd
import csv 

#topdir='/directory/where/experimental/data/is/located/'
topdir='./RigidLibrary/DataLattice'

# experimental friction coefficient
mu=5

# Experiment type, for other types please use a different runscript.
datatype = 'lattice'
bc='open'
bc='x' 
#bc='y'
#bc='xy'

nhex = 100 ########## please put the lattice size
# Potential loop over multiple lattices
lattice_nums=['_'+str(nhex)+'by'+str(nhex)+'_']

# Loop over experiment
for lattice in lattice_nums:
        # Extract some information from the data set:
        filename =topdir+'/Adjacency_list' + lattice + bc + '.txt'
        print(filename)
        try:
            data = np.loadtxt(topdir+'/Adjacency_list' + lattice + bc + '.txt', delimiter=',')
            #nsteps is the maximal frame number
            nsteps = np.max(data[:,0]).astype(int)
        except:
            print('No valid data found for lattice ' + lattice)
            #if no data is found for an experiment number go to the next experiment
            continue
        
        
        #Creating configuration
        ThisConf = CF.Configuration(topdir,datatype, bc, nhex)
                
        #Reading in the data
        ThisConf.readLatticedata(lattice, bc,True)

        # Play the pebble game
        ########## Setting up and playing the pebble game
        ThisPebble = PB.Pebbles(ThisConf,3,3,'nothing',False)

        # play game
        ThisPebble.play_game()
        # compute rigid clusters
        cidx, clusterall, clusterallBonds, clusteridx, BigCluster=ThisPebble.rigid_cluster()

        print(ThisPebble.pcluster)

        maxlabels = ThisPebble.MaxRigidPos()
        pdata = np.zeros((len(maxlabels),2))
        pdata[:,0] = ThisConf.x[maxlabels]
        pdata[:,1] = ThisConf.y[maxlabels]
        filename = topdir + '/MaxCluster_' + bc + '.dat'
        with open(filename,'w',newline="\n") as f:
            wr=csv.writer(f)
            wr.writerows(pdata)

        ########## Have a look at some analysis functions of the rigid clusters
        #def __init__(self,conf0,pebbles0,hessian0,tiling0='skip',fgiven0=0.001,verbose0=False):
        ThisAnalysis=AN.Analysis(ThisConf,ThisPebble,0)

        #def plotStresses(self,plotCir,plotVel,plotCon,plotF,plotStress,**kwargs):
        fig1 = ThisAnalysis.plotStresses(False,False,True,False,False)    

        #def plotPebbles(self,plotCir,plotPeb,plotPebCon,plotClus,plotOver,**kwargs):  
        fig2 = ThisAnalysis.plotPebbles(True,True,True,False,True)  

        fig3 = ThisAnalysis.plotPebbles(True,True,False,True,False)  
        
        plt.show()