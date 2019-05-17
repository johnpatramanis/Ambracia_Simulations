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
from multiprocessing import Process,Manager


start_time = time.time()


REPS=1


##########################################################################################################################################


N_Boet1i=100
N_Boet2i=100
N_Boet3i=100

N_Foreign1=10000
N_Foreign2=10000

generation_time = 20

T=500/generation_time
T_old=3000/generation_time
T_old_old=5000/generation_time



r_boet=0.1


N_Boet1= N_Boet1i/ math.exp(-r_boet * (T/generation_time))
N_Boet2= N_Boet2i/ math.exp(-r_boet * (T/generation_time))
N_Boet3= N_Boet3i/ math.exp(-r_boet * (T/generation_time))



###############################################################################################################################



population_configurations = [
    msprime.PopulationConfiguration(initial_size=N_Boet1i,growth_rate=r_boet),
    msprime.PopulationConfiguration(initial_size=N_Boet2i,growth_rate=r_boet),
    msprime.PopulationConfiguration(initial_size=N_Boet3i,growth_rate=r_boet),
    msprime.PopulationConfiguration(initial_size=N_Foreign1,growth_rate=0),
    msprime.PopulationConfiguration(initial_size=N_Foreign2,growth_rate=0)

]



migration_matrix = [
[0,0.001,0.001,0.001,0.001],
[0.001,0,0.001,0.01,0.01],
[0.001,0.001,0,0.1,0.1],
[0.001,0.01,0.1,0,0.001],
[0.001,0.01,0.1,0.001,0]
]

N1=50
N2=50
N3=50




POPS=[N1,N2,N3]
samples=[msprime.Sample(0,0)]*N1 + [msprime.Sample(1,0)]*N2 + [msprime.Sample(2,0)] *N3




demographic_events = [

msprime.PopulationParametersChange(time =T , growth_rate = 0 , population_id = 0),
msprime.PopulationParametersChange(time =T ,growth_rate=0 , population_id = 1),
msprime.PopulationParametersChange(time =T ,growth_rate=0 , population_id = 2),

msprime.PopulationParametersChange(time =T , initial_size = N_Boet1i , population_id = 0),
msprime.PopulationParametersChange(time =T , initial_size = N_Boet2i, population_id = 1),
msprime.PopulationParametersChange(time =T , initial_size = N_Boet3i, population_id = 2),

msprime.MigrationRateChange(time=T , rate=0, matrix_index=(0,3)), 
msprime.MigrationRateChange(time=T , rate=0, matrix_index=(3,0)),
msprime.MigrationRateChange(time=T , rate=0, matrix_index=(4,0)),
msprime.MigrationRateChange(time=T , rate=0, matrix_index=(0,4)),


msprime.MassMigration(time=T,source=2,destination=1,proportion = 1.0),
msprime.MassMigration(time=T,source=1,destination=0,proportion = 1.0),
msprime.MassMigration(time=T_old,source=3,destination=0,proportion = 1.0),
msprime.MassMigration(time=T_old_old,source=4,destination=0,proportion = 1.0)




]






##################################################################################################################################################


def SIMULATE(L,argument,samples,population_configurations,migration_matrix,demographic_events):
    j=int(argument)
    #recomb_map=msprime.RecombinationMap.read_hapmap('genetic_map_GRCh37_chr{}.txt'.format(j))
    dd = msprime.simulate(samples=samples,
        population_configurations=population_configurations,
        migration_matrix=migration_matrix,mutation_rate=1e-8,recombination_rate=2e-8,
        demographic_events=demographic_events,length=150000000)
    outfile=open('ms_prime_{}'.format(j),'w')
    for var in dd.variants():
        L.append([int(j),var.index,var.position])
        for genotype in var.genotypes:
            outfile.write(str(genotype))
        outfile.write('\n')
    outfile.close()    
    wow=open('mynewvcf{}.vcf'.format(j),'w')
    dd.write_vcf(wow,2,str(j))
    wow.close()
    
    population_labels= ["Boet1"]*int(N1/2) + ["Boet2"]*int(N2/2) + ["Boet3"]*int(N3/2)
    d=0
    newlabels=[]
    for i in range(0,len(population_labels)):
        newlabels.append(population_labels[i]+str(d))
        d+=1
        if i==len(population_labels)-2:
            newlabels.append(population_labels[i]+str(d))
            break
        if population_labels[i]!=population_labels[i+1]:
            d=0
    population_labels=newlabels
    wow=open('mynewvcf{}.vcf'.format(j))
    wowzers=open('myvcf{}.vcf'.format(j),'w')
    for line in wow:
        line=line.strip().split()
        if line[0]=='#CHROM':
            line[9:]=population_labels
        wowzers.write("\t".join(line))
        wowzers.write("\n")
    wow.close()
    
    return j,L
L=[]
SIMULATE(L,0,samples,population_configurations,migration_matrix,demographic_events)



elapsed_time_1 = time.time() - start_time        
    
print('Step 1 : {} '.format(elapsed_time_1/60))        
