import msprime
import numpy as np
import math
import os
import re



N_OG=1000
N_OUT=1000
N_AB=1000
N_A0=1000
N_B0=1000

r_A=0.000
r_B=0.000

generation_time = 25

T_split_OUT_AB=25000/generation_time
T_split_AB=12500/generation_time




N_A=N_A0 / math.exp(-r_A * T_split_AB)
N_B=N_B0 / math.exp(-r_B * T_split_AB)



population_configurations = [
    msprime.PopulationConfiguration(initial_size=N_OUT),
    msprime.PopulationConfiguration(initial_size=N_A, growth_rate=r_A),
    msprime.PopulationConfiguration(initial_size=N_B, growth_rate=r_B)
]

#migration of AB to OUT
m_OUT_AB=0.00000

migration_matrix = [
    [0,0.00000,0.000000],
    [0.00000,0,0.00000],
    [0.000000,0.00000,0]]

samples=[msprime.Sample(0,0)]*100 + [msprime.Sample(1,0)]*100 + [msprime.Sample(2,0)] *100


demographic_events = [
    # A and B merge
    msprime.MassMigration(time=T_split_AB, source=2, destination=1, proportion=1.0),
    msprime.MigrationRateChange(time=T_split_AB, rate=0),
    msprime.MigrationRateChange(time=T_split_AB, rate=m_OUT_AB, matrix_index=(0, 1)),
    msprime.MigrationRateChange(time=T_split_AB, rate=m_OUT_AB, matrix_index=(1, 0)),
    msprime.PopulationParametersChange(time=T_split_AB, initial_size=N_AB, growth_rate=0, population_id=1),
    # Population AB merges into OUT
    msprime.MassMigration(time=T_split_OUT_AB, source=1, destination=0, proportion=1.0),
    # Size changes to N_A at T_AF
    msprime.PopulationParametersChange(time=T_split_OUT_AB, initial_size=N_OUT, population_id=0)
]

dd = msprime.simulate(samples=samples,
    population_configurations=population_configurations,
    migration_matrix=migration_matrix,mutation_rate=1e-8,
    demographic_events=demographic_events,length=1000000, recombination_rate=2e-8,num_replicates=22)

    
MYRUN=1
for RUN in dd:
    variants=[]
    outfile=open('ms_prime_{}'.format(MYRUN),'w')
    for variant in RUN.variants():
        #print(variant)
        variants.append([variant.position,variant.index])
        for genotype in variant.genotypes:
            outfile.write(str(genotype))
        outfile.write('\n')
    MYRUN+=1
    outfile.close()
    
################################## 
#FLIP GENOTYPES
MAXRUNS=MYRUN
MYRUN=1
while MYRUN<MAXRUNS:
    
    msfile=open('ms_{}'.format(MYRUN),'w')
    column=0
    while column<len(samples):
        msprimefile=open('ms_prime_{}'.format(MYRUN),'r')
        person=[]
        for line in msprimefile:
            line=line.strip().split()[0]
            person.append(str(line[column]))
        msfile.write(''.join(person))
        msfile.write('\n')
        column+=1
        msprimefile.close()
    MYRUN+=1

#####################################################
#MERGING OF FILES ####
##
REPLICATE=0
MERGED=open('ms_allchroms_{}'.format(REPLICATE),'w')
for sample in range(0,len(samples)):
    for chromosome in range(1,23):
        file=open('ms_{}'.format(chromosome),'r')
        myline=0
        for line in file:
            #print(sample,myline)
            if myline==sample:
                line=line.strip()
                MERGED.write(str(line))
                break
            myline+=1
        file.close()
    MERGED.write('\n')
MERGED.close()

## os.system('rm ms_*.')
