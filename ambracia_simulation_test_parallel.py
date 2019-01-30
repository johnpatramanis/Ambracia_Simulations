import msprime
import numpy as np
import math
import os
import time
import re
import random
from multiprocessing import Process


start_time = time.time()


reps=1
for REPS in range(0,reps):

    totalf3=[]
    

    N_OG=1000
    N_OUT=1000
    N_AB=1000
    N_A0=1000
    N_B0=1000

    r_A=0.000
    r_B=0.000

    generation_time = 25

    T_split_OUT_AB=10000/generation_time
    T_split_AB=5000/generation_time




    N_A=N_A0 / math.exp(-r_A * T_split_AB)
    N_B=N_B0 / math.exp(-r_B * T_split_AB)



    population_configurations = [
        msprime.PopulationConfiguration(initial_size=N_OUT),
        msprime.PopulationConfiguration(initial_size=N_A, growth_rate=r_A),
        msprime.PopulationConfiguration(initial_size=N_B, growth_rate=r_B)
    ]

    #migration of AB to OUT
    m_OUT_AB=0.001

    migration_matrix = [
        [0,0.0001,0.0001],
        [0.0001,0,0.0001],
        [0.0001,0.0001,0]]

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

    
    
    
    
    




######################################################################################################################################################
#RUN the simulation and output genotypes in vcfs and ms format files, one for each chrom 
    variantinfo=[]
    variants=[]
    
    def SIMULATE(argument):
        j=int(argument)
        recomb_map=msprime.RecombinationMap.read_hapmap('genetic_map_GRCh37_chr{}.txt'.format(j))
        dd = msprime.simulate(samples=samples,
            population_configurations=population_configurations,
            migration_matrix=migration_matrix,mutation_rate=1e-8,
            demographic_events=demographic_events,recombination_map=recomb_map)
        outfile=open('ms_prime_{}'.format(j),'w')   
        for var in dd.variants():
            variants.append([var.index,var.position])
            variantinfo.append('{}\t{}\t{}\n'.format(j,var.index,var.position))
            for genotype in var.genotypes:
                outfile.write(str(genotype))
            outfile.write('\n')
        outfile.close()    
        wow=open('mynewvcf{}.vcf'.format(j),'w')
        dd.write_vcf(wow,2,str(j))
        wow.close()
        
        population_labels= ["africa"]*50 + ["asia"]*50 + ["europe"]*50
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
        
        return j
    
    
    
    
    if __name__ == '__main__':
        processes=[]
        for loop in range(1,23):
            p=Process(target=SIMULATE,args=(loop,))
            processes.append(p)
            p.start()
    
            
        for p in processes:
            p.join()
    

    
    variantinfo=sorted(variantinfo)
    variantinformation=open('variants_info.txt','w')
    variantinformation.write('CHROM\tVARIANT\tPOSITION\n')
    for loop in variantinfo:
        variantinformation.write(loop)
    
    variantinformation.close()

    elapsed_time_1 = time.time() - start_time        
        
    print('Step 1 : {} '.format(elapsed_time_1/60))        
        
        
        
        
        
######################################################################################################################################################
#Transform msprime format files to ms format
#prepare for COMUS stats

    MYRUN=22
    MAXRUNS=MYRUN
    MYRUN=1
    while MYRUN<=MAXRUNS:
        
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

######################################################

#MERGING OF FILES ####


    MERGED=open('ms_allchroms_{}'.format(REPS),'w')
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
    
#####################################################
#Split each ms format chromosome file to 50kb chunks

    SNPS=open('variants_info.txt','r')
    firstLine = SNPS.readline()
    CHUNKS=[]
    POSITIONS=[]
    begin=0
    end=0
    counter=0
    chr=1
    
    for line in SNPS:
        line=line.strip().split()
        POSITIONS.append(line[2])
        end=float(line[2])
        if (end-begin>=500000.0) :
            CHUNKS.append([counter,begin,end])
            begin=float(line[2])
        if chr!=int(line[0]):
            begin=0
            chr=int(line[0])
            
        counter+=1
    SNPS.close()
    MS_MERGED=open('ms_allchroms_{}'.format(REPS),'r')    
    

    
    MS_ALL_CHROMS=[]
    for line in MS_MERGED:        
        line=line.strip().split()
        MS_ALL_CHROMS.append(line[0])
    #print(len(MS_ALL_CHROMS))
    #print(type(MS_ALL_CHROMS))
    counter=0
    begin=0
    opener=open('CHUNKED_{}'.format(REPS),'w')
    opener.write('ms {} {}\n{} {} {}'.format(len(samples),len(CHUNKS),random.randint(0,10000),random.randint(0,10000),random.randint(0,10000)))
    opener.write('\n')
    for x in CHUNKS:
        opener.write("\n")
        opener.write("//\n")
        opener.write("segsites: {}\n".format(x[0]-begin))
        positions_of_this_chunk=POSITIONS[begin:x[0]]
        positions_of_this_chunk=' '.join(positions_of_this_chunk)
        opener.write("positions: {}".format(positions_of_this_chunk))
        opener.write("\n")

        for y in MS_ALL_CHROMS:                    
            opener.write(''.join(y[begin:x[0]]))
            opener.write('\n')            
        begin=x[0]+1
        
        counter+=1

        
    opener.close()
    elapsed_time_2=time.time() - start_time
    print('step 2 : {}'.format(elapsed_time_2/60))



#Arrange the vcf files into one, fix labels , bed file , transform to eigen, calculate pca and f stats
        
    os.system('rm mynewvcf*.vcf')
    os.system('bcftools concat -o total_chroms.vcf myvcf*.vcf')
    os.system('rm myvcf*.vcf')




    VCF=open('total_chroms.vcf','r')
    newVCF=open('newtotal_chroms.vcf','w')

    snpcount=0
    
    variants=sorted(variants)
    for line in VCF:
        if line[0]!='#' and snpcount<len(variants):
            line=line.strip().split()
            line[2]='rs{}'.format(snpcount)
            line[1]=str(variants[snpcount][1])
            line.append('\n')
            line='\t'.join(line)
            snpcount+=1
        newVCF.write(line)

    VCF.close
    newVCF.close

    os.system('mv newtotal_chroms.vcf total_chroms.vcf')

    os.system('plink --vcf total_chroms.vcf --make-bed --out simulation')



    

    import os.path
    if os.path.isfile('simulation.bed'):
        simulationfile='simulation'
    else:
        simulationfile='simulation-temporary'
    
    os.system('plink --bfile {} --pca 10 --out pcaofsimulation'.format(simulationfile))
    
    
    ####################################### 3 Pop Test ######################################################################################
    parfile=open('parfile.txt','w')

    parfile.write('genotypename:    {}.bed\n'.format(simulationfile))
    parfile.write('snpname:         {}.bim\n'.format(simulationfile))
    parfile.write('indivname:       {}.fam\n'.format(simulationfile))
    parfile.write('outputformat:   PACKEDANCESTRYMAP\n')
    parfile.write('genotypeoutname: simulation.geno\n')
    parfile.write('snpoutname:      simulation.snp\n')
    parfile.write('indivoutname:    simulation.ind\n')
    parfile.write('pordercheck: NO')

    parfile.close()


    os.system('convertf -p parfile.txt')

    IND=open('simulation.ind','r')
    newIND=open('newsimulation.ind','w')

    for line in IND:
        line=line.strip().split()
        label=re.search(r'([a-z]+)([0-9]+):[a-z]+[0-9]+',line[0])
        pop=label.group(1)
        number=label.group(0)
        line[0]=str(number)
        line[2]=str(pop)
        newIND.write('\t'.join(line))
        newIND.write('\n')
        
    IND.close
    newIND.close()
        
    os.system('mv newsimulation.ind simulation.ind')

    Pop3=open('qp3Poplist','w')
    Pop3.write('africa europe asia')
    Pop3.close()

    Parfilepop=open('3popparfile','w')
    Parfilepop.write('SSS: allmap\n')
    Parfilepop.write('indivname:   simulation.ind\n')
    Parfilepop.write('snpname:     simulation.snp\n')
    Parfilepop.write('genotypename: simulation.geno\n')
    Parfilepop.write('popfilename: qp3Poplist\n')

    Parfilepop.close()


    #SNP=open('simulation.snp','r')
    #newSNP=open('newsimulation.snp','w')
    #snpcounter=0
    #for line in SNP:
    #    if snpcounter<len(variants):
    #        line=line.strip().split()
    #        line[0]='rs{}'.format(snpcounter)
    #        line[2]=str(variants[snpcounter])
    #        line.append('\n')
    #        line='\t'.join(line)
    #        snpcounter+=1
    #        newSNP.write(line)




    #SNP.close
    #newSNP.close

    os.system('mv newsimulation.snp simulation.snp')
    os.system('ls')



    os.system('qp3Pop -p 3popparfile >f3stat_{}'.format(REPS))
    
    
    f3file=open('f3stat_{}'.format(REPS),'r')
    
    for line in f3file:
        line=line.strip().split()
        #print(line)
        if line[0]=='result:':
            totalf3.append(float(line[4]))
    f3file.close()
    
    
    
    
    os.system('rm simulation.*')
    os.system('rm simulation-temporary.*')
    os.system('rm ms_prime_*')
    for x in range(1,23):
        os.system('rm ms_{}'.format(x))

    elapsed_time_3=time.time() - start_time
    print('step 3 : {}'.format(elapsed_time_3/60))        
    
    os.system('CoMuStats -input CHUNKED_{} -npop 3 100 100 100 -ms > COMUSTATS_{}'.format(REPS,REPS))
    elapsed_time_4=time.time() - start_time
    print('step 4 : {}'.format(elapsed_time_4/60)) 
        
###############################################################################################################################################
  
    

###############################################################################################################################################
print(np.mean(totalf3))

