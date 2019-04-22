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


N_Boet1=10000
N_Boet1_High_soc=500

N_Boet2=5000
N_Boet3=5000
N_Boet4=5000

N_Foreign1=10000
N_Foreign2=10000

generation_time = 20

r_boet=0.0001



###############################################################################################################################



population_configurations = [
    msprime.PopulationConfiguration(initial_size=N_Boet1,growth_rate=r_boet),
    msprime.PopulationConfiguration(initial_size=N_Boet1_High_soc,growth_rate=r_boet),
    msprime.PopulationConfiguration(initial_size=N_Boet2,growth_rate=r_boet),
    msprime.PopulationConfiguration(initial_size=N_Boet3,growth_rate=r_boet),
    msprime.PopulationConfiguration(initial_size=N_Boet4,growth_rate=r_boet),
    msprime.PopulationConfiguration(initial_size=N_Foreign1,growth_rate=r_boet),
    msprime.PopulationConfiguration(initial_size=N_Foreign2,growth_rate=r_boet)

]



migration_matrix = [
    [0,0.001,0.0001,0.0001,0.0001,0.00001,0.00001],
    [0.001,0.0,0.0,0.0,0.0,0.0,0.0],
    [0.0001,0.0,0.0,0.0001,0.0001,0.00001,0.00001],
    [0.0001,0.0,0.0001,0.0,0.0001,0.00001,0.00001],
    [0.0001,0.0,0.0001,0.0001,0.0,0.00001,0.00001],
    [0.00001,0.0,0.00001,0.00001,0.00001,0.0,0.00001],
    [0.00001,0.0,0.00001,0.00001,0.00001,0.00001,0.0]

]

N1=50
N2=20
N3=50
N4=50
N5=50
N6=20
N7=20

POPS=[N1,N2,N3,N4,N5,N6,N7]
samples=[msprime.Sample(0,0)]*N1 + [msprime.Sample(1,0)]*N2 + [msprime.Sample(2,0)] *N3 +[msprime.Sample(3,0)]*N4 + [msprime.Sample(4,0)]*N5 + [msprime.Sample(5,0)] *N6 + [msprime.Sample(6,0)] *N7




demographic_events = [

]






##################################################################################################################################################


def SIMULATE(L,argument,samples,population_configurations,migration_matrix,demographic_events):
    j=int(argument)
    #recomb_map=msprime.RecombinationMap.read_hapmap('genetic_map_GRCh37_chr{}.txt'.format(j))
    dd = msprime.simulate(samples=samples,
        population_configurations=population_configurations,
        migration_matrix=migration_matrix,mutation_rate=1e-8,
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
    
    population_labels= ["Boet1"]*int(N1/2) + ["Boet1_High_Soc"]*int(N2/2) + ["Boet2"]*int(N3/2) +["Boet3"]*int(N4/2) + ["Boet4"]*int(N5/2) + ["Foreign1"]*int(N6/2) +  ["Foreign2"]*int(N7/2)
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
