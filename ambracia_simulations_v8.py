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


reps=1000
for REPS in range(0,reps):


##############################################################################################################################################
#Simulation Parameters
    
    parametersfile=open('PARAMETERS_{}'.format(REPS),'w')
    
    N_locals=int(round(random.uniform(500.0,1000.0)))
    N_metropolis=int(round(random.uniform(500.0,1000.0)))
    N_outgroup=int(round(random.uniform(800.0,1200.0)))
    
    generation_time = 20
    T_COLONIZATION=700/generation_time
    
    
    COLONIZER=random.randint(0,1)
    if COLONIZER==0:
        N_initial_colony=int(round(random.uniform(200.0,float(N_locals))))
        while N_initial_colony>N_metropolis:
            N_initial_colony=int(round(random.uniform(200.0,float(N_metropolis))))
    if COLONIZER==1:
        N_initial_colony=int(round(random.uniform(200.0,float(N_metropolis))))



    r_locals=10**(-1*random.uniform(1,4))
    r_metropolis=10**(-1*random.uniform(1,4))
    r_colony=10**(-1*random.uniform(1,4))
    
    print(N_locals,N_metropolis,N_initial_colony,r_locals,r_metropolis)
    
    growth_counter=0
    while (float(N_initial_colony) / (math.exp(-r_colony * T_COLONIZATION)) ) > float(N_metropolis):
        r_colony=10**(-1*random.uniform(1,4))
        growth_counter+=1
        if growth_counter>=1000000:
            r_colony=0
            break
            
    N_finale_colony=N_initial_colony / (math.exp(-r_colony * T_COLONIZATION))
    
##################################################################################################################################
#### MSPRIME SET-UP


    population_configurations = [
        msprime.PopulationConfiguration(initial_size=N_locals,growth_rate=r_locals),
        msprime.PopulationConfiguration(initial_size=N_metropolis, growth_rate=r_metropolis),
        msprime.PopulationConfiguration(initial_size=N_finale_colony, growth_rate=r_colony),
        msprime.PopulationConfiguration(initial_size=N_outgroup, growth_rate=0)
    ]



    migration_matrix = [
        [0,10**(-1*random.uniform(1,4)),10**(-1*random.uniform(1,4)),0.0001],
        [10**(-1*random.uniform(1,4)),0,10**(-1*random.uniform(1,4)),0.0001],
        [10**(-1*random.uniform(1,4)),10**(-1*random.uniform(1,4)),0,0.0001],
        [0.0001,0.0001,0.0001,0]
        ]

    N1=20
    N2=20
    N3=20
    N4=20
    POPS=[N1,N2,N3]
    samples=[msprime.Sample(0,0)]*N1 + [msprime.Sample(1,0)]*N2 + [msprime.Sample(2,0)] *N3 + [msprime.Sample(3,0)]*N4

    demographic_events = [
    msprime.MigrationRateChange(time=T_COLONIZATION, rate=0, matrix_index=(0, 2)),
    msprime.MigrationRateChange(time=T_COLONIZATION, rate=0, matrix_index=(2, 0)),
    msprime.MigrationRateChange(time=T_COLONIZATION, rate=0, matrix_index=(1, 2)),
    msprime.MigrationRateChange(time=T_COLONIZATION, rate=0, matrix_index=(2, 1)),
    msprime.MigrationRateChange(time=T_COLONIZATION, rate=0, matrix_index=(2, 3)),
    msprime.MigrationRateChange(time=T_COLONIZATION, rate=0, matrix_index=(3, 2)),
    
    msprime.PopulationParametersChange(time=T_COLONIZATION, initial_size=N_initial_colony, growth_rate=0, population_id=2),
    msprime.MassMigration(time=T_COLONIZATION, source=2, destination=COLONIZER, proportion=1.0),
    
    
    ]

    parametersfile.write('\t'.join([str(x) for x in [COLONIZER,N_locals,N_metropolis,N_initial_colony,N_finale_colony,r_locals,r_metropolis,r_colony]]))
    parametersfile.write('\n')
    parametersfile.write('\t'.join([str(x) for x in migration_matrix]))
    
    print(migration_matrix)
    
    




######################################################################################################################################################
#RUN the simulation and output genotypes in vcfs and ms format files, one for each chrom 

    
    def SIMULATE(L,argument,samples,population_configurations,migration_matrix,demographic_events):
        j=int(argument)
        recomb_map=msprime.RecombinationMap.read_hapmap('genetic_map_GRCh37_chr{}.txt'.format(j))
        dd = msprime.simulate(samples=samples,
            population_configurations=population_configurations,
            migration_matrix=migration_matrix,mutation_rate=1e-8,
            demographic_events=demographic_events,recombination_map=recomb_map)
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
        
        population_labels= ["locals"]*int(N1/2) + ["metropolis"]*int(N2/2) + ["apoikia"]*int(N3/2) + ["outgroup"]*int(N3/2)
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
    
    
###############################################################################################################################################
# MULTITHREDING #
    L=[]
    if __name__ == '__main__':
        with Manager() as manager:
            L=manager.list(L)
            processes=[]
            for loop in range(1,23):
                p=Process(target=SIMULATE,args=(L,loop,samples,population_configurations,migration_matrix,demographic_events,))
                processes.append(p)
                
                p.start()
        
                
            for p in processes:
                p.join()
            #print(len(L),'1')
            sys.stdout.flush()
            variants=sorted(list(L))


    variantinfo=['{}\t{}\t{}\n'.format(x[0],x[1],x[2])for x in variants]
    print(len(variants),len(variantinfo))

    variantinformation=open('variants_info.txt','w')
    variantinformation.write('CHROM\tVARIANT\tPOSITION\n')
    for loop in variantinfo:
        variantinformation.write(loop)
    
    variantinformation.close()

    elapsed_time_1 = time.time() - start_time        
        
    print('Step 1 : {} '.format(elapsed_time_1/60))        
    print(REPS)
        
        
        
        
######################################################################################################################################################
#Transform msprime format files to ms format
#prepare for COMUS stats

    MYRUN=22
    MAXRUNS=MYRUN
    MYRUN=1
    while MYRUN<=MAXRUNS:
        
        msfile=open('ms_{}'.format(MYRUN),'w')
        column=0
        while column<len(samples)-N4:
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

##############################################################################################################################################

#MERGING OF FILES ####


    MERGED=open('ms_allchroms_{}'.format(REPS),'w')
    for sample in range(0,len(samples)-N4):
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
    
#############################################################################################################################################
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
########################################## CHUNKS TO BE REMOVED? #############################################################################
    TO_BE_REMOVED=[]
    
    for x in range(0,len(CHUNKS)-1):
        #if ((CHUNKS[x+1][1] - CHUNKS[x][2]) <= 500000.0) and (CHUNKS[x] not in TO_BE_REMOVED) and (CHUNKS[x+1][1] - CHUNKS[x][2] >= 0):
        if x%2!=0:
            TO_BE_REMOVED.append((CHUNKS[x+1]))
    
    
    
    
    
    
    
    
    
    
    MS_MERGED=open('ms_allchroms_{}'.format(REPS),'r')    
#################################### PRINT CHUNKS MS FORMAT FOR COMUSTATS ####################################################################
    MS_ALL_CHROMS=[]
    for line in MS_MERGED:        
        line=line.strip().split()
        MS_ALL_CHROMS.append(line[0])
    counter=0
    begin=0
    opener=open('CHUNKED_{}'.format(REPS),'w')
    opener.write('ms {} {}\n{} {} {}'.format((len(samples)-N4),len(CHUNKS),random.randint(0,10000),random.randint(0,10000),random.randint(0,10000)))
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
##################################REPRINT WITH GAPS OF 500kb BETWEEN CHUNKS ######################################################################
        
    opener.close()
    elapsed_time_2=time.time() - start_time
    opener=open('CHUNKED_{}'.format(REPS),'r')
    opener2=open('CHUNKED_500kb_gaps_{}'.format(REPS),'w')
    line_counter=0
    for line in opener:
        if line=="//\n":
            line_counter+=1
        line=line.split()
        if line_counter%2==0:
            opener2.write(' '.join(line))
            opener2.write('\n')
    
    
    opener.close()
    opener2.close()
    print('step 2 : {}'.format(elapsed_time_2/60))



#############Arrange the vcf files into one, fix labels , bed file , transform to eigen, calculate pca and f stats ############################
        
    os.system('rm mynewvcf*.vcf')
    os.system('bcftools concat -o total_chroms.vcf myvcf*.vcf')
    os.system('rm myvcf*.vcf')




    VCF=open('total_chroms.vcf','r')
    newVCF=open('newtotal_chroms.vcf','w')
    VCFinfo={}
    snpcount=0
    
    variants=sorted(variants)
    for line in VCF:
        if line[0]!='#' and snpcount<len(variants):
            line=line.strip().split()
            if len(line)<=2:
                continue
            line[2]='rs{}'.format(snpcount)
            line[1]=str(variants[snpcount][2])
            if line[0] in VCFinfo:
                VCFinfo[line[0]].append([float(line[1]),line[2]])
            else:
                VCFinfo[line[0]]=[[float(line[1]),line[2]]]
            line.append('\n')
            line='\t'.join(line)
            snpcount+=1
        newVCF.write(line)

    VCF.close
    newVCF.close

    os.system('mv newtotal_chroms.vcf total_chroms.vcf')

    os.system('plink --vcf total_chroms.vcf --make-bed --out simulation')


    

    if os.path.isfile('simulation.bed'):
        simulationfile='simulation'
    else:
        simulationfile='simulation-temporary'
    
    
    removefam=open('famtoremove.txt','w')
    for k in range(0,int((N4/2)-1)):
        removefam.write('outgroup{}\n'.format(k))
    removefam.close()
    
    os.system('plink --bfile {} --remove-fam famtoremove.txt --pca 10 --out pcaofsimulation'.format(simulationfile))
    
    
    PCAFILE=open('pcaofsimulation.eigenvec','r')
    eigenvecs=[]
    for line in PCAFILE:
        line=line.strip().split()
        eigenvecs.append([float(line[2]),float(line[3])])

    counter=0
    pca_clust=[]    
    for pop in POPS:
        stop=int(pop/2)+counter
        print(counter)
        print(stop)
        pca_clust.append(eigenvecs[counter:stop])
        counter=stop

    print(pca_clust)

    indistance=[]
    centers=[]
    outdistance=[]

    for cluster in pca_clust:
        clustercenter=np.mean(np.asarray(cluster),axis=0)
        centers.append(clustercenter)
        meandistfromcent=[]
        for x in cluster:
            meandistfromcent.append(numpy.linalg.norm(clustercenter-x))
        meandistfromcent=np.mean(meandistfromcent)
        indistance.append(meandistfromcent)

    for j in centers:
        for k in centers:
            if numpy.linalg.norm(j-k)!=0 and numpy.linalg.norm(j-k) not in outdistance:
                outdistance.append(numpy.linalg.norm(j-k))

    print(indistance)
    print(outdistance)

    PCACLUSTERING=open('PCA_CLUSTERING_{}'.format(REPS),'w')
    PCACLUSTERING.write('\t'.join([str(x) for x in indistance])+'\t'+'\t'.join([str(x) for x in outdistance]))
    PCACLUSTERING.close()


    
    
    
    
    
    
    
############################################# 3 Pop Test ######################################################################################
    for k in range(1,22):
        begin=0
        end=0
        segments=[]
        segments.append(VCFinfo[str(k)][0][1])
        for y in VCFinfo[str(k)]:
            end=float(y[0])
            if ((end-begin)>=10000000.0):
                segments.append(y[1])
                begin=float(y[0])+100000.0
    totalf3=[]
    for j in range(0,len(segments)-1):
        print(segments[j])
        os.system('plink --vcf total_chroms.vcf  --from {} --to {} --make-bed --out simulation'.format(segments[j],segments[j+1]))

        if os.path.isfile('simulation.bed'):
            simulationfile='simulation'
        else:
            simulationfile='simulation-temporary'
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
        SNP=open('simulation.snp','r')
        snpcounter=0
        for line in SNP:
            snpcounter+=1
        SNP.close()
        
        
        os.system('mv newsimulation.ind simulation.ind')
    
        Pop3=open('qp3Poplist','w')
        Pop3.write('locals metropolis apoikia')
        Pop3.close()
    
        Parfilepop=open('3popparfile','w')
        Parfilepop.write('SSS: allmap\n')
        Parfilepop.write('indivname:   simulation.ind\n')
        Parfilepop.write('snpname:     simulation.snp\n')
        Parfilepop.write('genotypename: simulation.geno\n')
        Parfilepop.write('popfilename: qp3Poplist\n')
    
        Parfilepop.close()
        os.system('mv newsimulation.snp simulation.snp')
        os.system('ls')
    
    
    
        os.system('qp3Pop -p 3popparfile >f3stat_{}'.format(REPS))
        
        
        f3file=open('f3stat_{}'.format(REPS),'r')
        for line in f3file:
            line=line.strip().split()
            #print(line)
            if line[0]=='result:':
                totalf3.append([float(line[4]),snpcounter])
        f3file.close()
        os.system('rm f3stat_{}'.format(REPS))
############################################## FINAL WRITING #################################################################################
    
    f3FINAL=open('f3FINAL_{}.txt'.format(REPS),'w')
    for line in totalf3:
        for x in line:
            f3FINAL.write(str(x))
            f3FINAL.write('\t')
        f3FINAL.write('\n')
    f3FINAL.close()
    
    os.system('rm famtoremove.txt')
    os.system('rm simulation.*')
    os.system('rm simulation-temporary.*')
    os.system('rm ms_prime_*')
    for x in range(1,23):
        os.system('rm ms_{}'.format(x))

    elapsed_time_3=time.time() - start_time
    print('step 3 : {}'.format(elapsed_time_3/60))        
    
    os.system('CoMuStats -input CHUNKED_500kb_gaps_{} -npop 3 20 20 20 -ms > COMUSTATS_{}'.format(REPS,REPS))
    elapsed_time_4=time.time() - start_time
    print('step 4 : {}'.format(elapsed_time_4/60)) 
        
###############################################################################################################################################
  
    

###############################################################################################################################################


