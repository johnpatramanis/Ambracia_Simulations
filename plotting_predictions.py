import numpy as np
import matplotlib.pyplot as plt


rep=1
file=open('PREDICTIONS_{}'.format(rep),'r')
actual=[]
predicted=[]
for line in file:
    line=line.strip().split()
    actual.append(line[0])
    predicted.append(line[1])


font = {'family': 'serif',
        'color':  'blue',
        'weight': 'normal',
        'size': 14,
        }
print(len(actual),len(predicted))
plt.figure()
plt.xlabel('Prediction', fontdict=font)    
plt.ylabel('Actual', fontdict=font)    
plt.title('Parameter number {}'.format(rep), fontdict=font)    
plt.plot(predicted,actual,'.', line)  
plt.savefig('prediction_{}.png'.format(rep))