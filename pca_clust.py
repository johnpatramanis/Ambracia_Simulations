import numpy as np





PCAFILE=open('ancient_ambrac.eigenvec','r')
eigenvecs=[]
for line in PCAFILE:
    line=line.strip().split()
    eigenvecs.append([float(line[2]),float(line[3])])


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
    counter=int(pop/2)

print(pca_clust)
for cluster in pca_clust:
    clustercenter=np.mean(np.asarray(cluster),axis=0)
    print(clustercenter)