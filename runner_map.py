import numpy as np
import sys
from scipy import integrate

sims=['BR','LCDM']

startPointNum = 3 

simulationTypes = 3


pairings=[]

for i in range(len(sims)):
    
    for j in range(startPointNum):
    
        for k in [2]:#range(simulationTypes):
    
            pairings.append([sims[i],j,k])

        

computerStreamNums=[4, 2]
streamNum=np.sum(computerStreamNums)
#streamNum must add up to len(sims)*startPointNum


index=0


scriptFiles=[]


for i in range(len(computerStreamNums)):

    runnerFile=open('./scripts3/runner_'+str(i)+'.sh', 'w')

    for j in range(computerStreamNums[i]):
        
        runnerFile.write('python3 ISW_map.py ' + pairings[index][0] + ' ' + str(pairings[index][1])+ ' ' + str(pairings[index][2]) + ' & \n')
        
        index+=1
        
    runnerFile.close()
        