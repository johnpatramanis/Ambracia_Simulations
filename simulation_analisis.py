import msprime
import numpy as np
import numpy.linalg
import math
import os
import time
import re
import random
import sys
import os.path


colonizers=[]
N_locals=[]
N_metropolis=[]
N_initial_colony=[]
N_final_colony=[]
r_locals=[]
r_metropolis=[]
r_colony=[]
PCA_in_dist_locals=[]
PCA_in_dist_metropolis=[]
PCA_in_dist_colony=[]
PCA_in_out_locals_metro=[]
PCA_in_out_locals_colony=[]
PCA_in_out_metropolis_colony=[]
MEAN_F3=[]
    
    
    
    
    

for REPS in range(0,500):
    ######################################################################################################################################3###
    #open files
    
    PARAM_FILE=open('PARAMETERS_{}'.format(REPS),'w')
    PCA_CLUST_FILE=open('PCA_CLUSTERING_{}'.format(REPS),'w')
    F3FILE=open('f3FINAL_{}'.format(REPS),'w')
    COMUS_FILE=open('COMUSTATS_{}'.format(REPS),'w')
    
    ###########################################################################################################################################
    #load parameters
    count=0
    for line in PARAM_FILE:
        line=line.strip().split()
        if count=0:
            colonizers.append(int(line[0]))
            N_locals.append(int(line[1]))
            N_metropolis.append(int(line[2]))
            N_initial_colony.append(int(line[3]))
            N_final_colony.append(float(line[4]))
            r_locals.append(float(line[5]))
            r_metropolis.append(float(line[6]))
            r_colony.append(float(line[7]))
            count+=1
        if count>0:
            pass
    
    for line in PCA_CLUST_FILE:
        line.strip().split()
        PCA_in_dist_locals.append(line[0])
        PCA_in_dist_metropolis.append(line[1])
        PCA_in_dist_colony.append(line[2])
        PCA_in_out_locals.append(line[3])
        PCA_in_out_locals_colony.append(line[4])
        PCA_in_out_metropolis_colony.append(line[5])
        
        
    
    f3mean=[]
    for line in F3FILE:
        line.strip().split()
        f3mean.append(float(line[0]))
    MEAN_F3.append(np.mean(f3mean))
    
    #for line in COMUS_FILE:
    #    line.strip().split()
    
    
    PARAM_FILE.close()
    PCA_CLUST_FILE.close()
    F3FILE.close()
    COMUS_FILE.close()
    
    
    
    
    
    
    
print(len(MEAN_F3),len(PCA_in_dist_colony),len(N_final_colony))
