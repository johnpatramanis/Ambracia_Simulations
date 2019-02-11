import numpy as np
import numpy.linalg




PCAFILE=open('pcaofsimulation.eigenvec','r')
eigenvecs=[]
for line in PCAFILE:
    line=line.strip().split()
    eigenvecs.append([float(line[2]),float(line[3])])

REPS=1
Pop1=20
Pop2=20
Pop3=20
POPS=[Pop1,Pop2,Pop3]
    
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
