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
from pydoc import help
from scipy.stats.stats import pearsonr




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
PCA_out_locals_metro=[]
PCA_out_locals_colony=[]
PCA_out_metropolis_colony=[]
MEAN_F3=[]
migration_local_colony=[]
migration_metropolis_colony=[]
migration_locals_metropolis=[]
migration_metropolis_locals=[]
    
    
    
    
    
    
    

for REPS in range(10,999):
    ######################################################################################################################################3###
    #open files
    
    PARAM_FILE=open('PARAMETERS_{}'.format(REPS),'r')
    PCA_CLUST_FILE=open('PCA_CLUSTERING_{}'.format(REPS),'r')
    F3FILE=open('f3stat_{}'.format(REPS),'r')
    COMUS_FILE=open('COMUSTATS_{}'.format(REPS),'r')
    
    ###########################################################################################################################################
    #load parameters
    
    for line in PARAM_FILE:
        if line[0]==str(1) or line[0]==str(0):
            line=line.strip().split()
            colonizers.append(int(line[0]))
            N_locals.append(int(line[1]))
            N_metropolis.append(int(line[2]))
            N_initial_colony.append(int(line[3]))
            N_final_colony.append(float(line[4]))
            r_locals.append(float(line[5]))
            r_metropolis.append(float(line[6]))
            r_colony.append(float(line[7]))
            
        else:
            line=line.strip()
            mysearch=re.findall(r'[0-9]\.[0-9]+',line)
            migration_locals_metropolis.append(float(mysearch[0]))
            migration_metropolis_locals.append(float(mysearch[2]))
            migration_local_colony.append(float(mysearch[4]))
            migration_metropolis_colony.append(float(mysearch[5]))
            
            #pop=migrations_numbers.group(1)
            #print(pop)
            
    
    for line in PCA_CLUST_FILE:
        line=line.strip().split('\t')
        PCA_in_dist_locals.append(float(line[0]))
        PCA_in_dist_metropolis.append(float(line[1]))
        PCA_in_dist_colony.append(float(line[2]))
        PCA_out_locals_metro.append(float(line[3]))
        PCA_out_locals_colony.append(float(line[4]))
        PCA_out_metropolis_colony.append(float(line[5]))
        
        
    
    f3mean=[]
    for line in F3FILE:
        line=line.strip().split()
            #print(line)
        if line[0]=='result:':
            f3mean.append(float(line[4]))
        
        
        
        
    MEAN_F3.append(np.mean(f3mean))
    
    #for line in COMUS_FILE:
    #    line.strip().split()
    
    
    PARAM_FILE.close()
    PCA_CLUST_FILE.close()
    F3FILE.close()
    COMUS_FILE.close()
    
    
    
    
parameters=[[migration_locals_metropolis,'migration_locals_from_metropolis'],[migration_metropolis_locals,'migration_metropolis_from_locals'],[migration_local_colony,'migration_local_from_colony'],[migration_metropolis_colony,'migration_metropolis_from_colony'],[colonizers,'who_colonized'],[N_locals,'N-locals'],[N_metropolis,'N_metropolis'],[N_initial_colony,'N_initial_colony'],[N_final_colony,'N_final_colony'],[r_locals,'r_locals'],[r_metropolis,'r_locals'],[r_colony,'r_colony']]
metrics=[[PCA_in_dist_locals,'PCA_in_dist_locals'],[PCA_in_dist_metropolis,'PCA_in_dist_metropolis'],[PCA_in_dist_colony,'PCA_in_dist_colony'],[PCA_out_locals_metro,'PCA_out_locals_metro'],[PCA_out_locals_colony,'PCA_out_locals_colony'],[PCA_out_metropolis_colony,'PCA_out_metropolis_colony'],[MEAN_F3,'F3_mean']]

CORRELATION=[]

for x in parameters:
    for y in metrics:
        CORRELATION.append([pearsonr(x[0],y[0]),x[1],y[1]])
        
CORELATION_FILE=open('CORRELATION_FILE','w')
CORELATION_FILE.write('Pearson_value'+'\t'+'P-value'+'\t'+'Parameter'+'\t'+'Metric'+'\n')


for writer in sorted(CORRELATION,key=lambda x: x[0][1]):
    CORELATION_FILE.write(str(writer[0][0])+'\t'+str(writer[0][1])+'\t'+str(writer[1])+'\t'+str(writer[2])+'\n')
