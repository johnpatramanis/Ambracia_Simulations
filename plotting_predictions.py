import numpy as np
import matplotlib.pyplot as plt
import math 

rep=1
file=open('PREDICTIONS_{}'.format(rep),'r')
actual=[]
predicted=[]
for line in file:
    line=line.strip().split()
    actual.append(float(line[0]))
    predicted.append(float(line[1]))

print(len(actual),len(predicted))
RMSE=[]
MAPE=[]

#Root Mean Squared Error 
for k in range(0,len(actual)):
    RMSE.append((actual[k]-predicted[k])**2)

ROOT=math.sqrt((sum(RMSE))/len(actual))
RMSE=ROOT/(max(actual)-min(actual))
print('Root Mean Squared Error is {} '.format(RMSE))

#Mean Absolute Percentage Error

for k in range(0,len(actual)):
    top=(actual[k]-predited[k])
    bot=(actual[k]+predited[k])/2
    
    MAPE.append(abs(top/bot) )
MAPE=sum(MAPE)/len(MAPE)

print('Mean Absolute Percentage Error is {} '.format(MAPE))



font = {'family': 'serif',
        'color':  'blue',
        'weight': 'normal',
        'size': 14,
        }

plt.figure()
plt.xlabel('Prediction', fontdict=font)    
plt.ylabel('Actual', fontdict=font)    
plt.title('Parameter number {}'.format(rep), fontdict=font)    
plt.plot(predicted,actual,'.', line)  
plt.savefig('prediction_{}.png'.format(rep))