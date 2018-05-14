import numpy as np
import sys
from scipy.interpolate import interp1d

startIndex=14
endIndex=64

computerStreamNums=[8,4,4]
streamNum=np.sum(computerStreamNums)

scriptFiles=[]

for i in range(len(computerStreamNums)):

    runnerFile=open('./scripts2/XXL_runner_'+str(i)+'.sh', 'w')

    for j in range(computerStreamNums[i]):
        scriptFiles.append(open('./scripts2/'+str(i)+'/XXL_'+str(j)+'.sh', 'w'))
        
        runnerFile.write('sh ./scripts2/'+str(i)+'/XXL_'+str(j)+'.sh &\n')
        
    runnerFile.close()


sims=['BR','LCDM']
 

exec(open('./loadCosmologies.py').read())

    
path1 = '../MilleniumXXL/'


outPath= './LC_results_XXL_lin/'


milleniumFiles=np.loadtxt(path1+'zlist_da.txt')


index=0

for i in range(len(sims)):
    
    tau_a_func=interp1d(scaleFactors[i], taus[i], kind='cubic', fill_value='extrapolate')
    
    for j in range(startIndex,endIndex):
        
        #WARNING: input snapshot has been fixed to the first one within AvERA cosmological curve coverage
        inFileName = path1+'density_1024_0'+str(14)+'.dat'      
        
        #scaleFactor=milleniumFiles[j,1]
        
        if (j>startIndex and j<endIndex-1):
        
            tauStart=str((tau_a_func(milleniumFiles[j-1,1])+tau_a_func(milleniumFiles[j,1]))/2.0)
            tauEnd=str((tau_a_func(milleniumFiles[j+1,1])+tau_a_func(milleniumFiles[j,1]))/2.0)
        
        elif (j==startIndex):
        
            tauStart=str(tau_a_func(milleniumFiles[j,1]))
            tauEnd=str((tau_a_func(milleniumFiles[j+1,1])+tau_a_func(milleniumFiles[j,1]))/2.0)
        
        elif (j==endIndex-1):
        
            tauStart=str((tau_a_func(milleniumFiles[j-1,1])+tau_a_func(milleniumFiles[j,1]))/2.0)
            tauEnd=str(tau_a_func(milleniumFiles[j,1]))

        tauInitial=str(tau_a_func(milleniumFiles[startIndex,1]))
        tauFinal=str(tau_a_func(milleniumFiles[endIndex-1,1]))
    
        outFileName = outPath+sims[i]+'/density_1024_0'+str(j)+'_'
        
        scriptFiles[index].write('python3 lightCurve_MXXL_lin.py ' + inFileName + ' ' + tauStart + ' ' + tauEnd + ' ' + tauInitial + ' ' + tauFinal + ' ' + outFileName + ' ' + sims[i] + ' \n')
        
        index=(index+1)%streamNum
        
        
for i in range(len(scriptFiles)):

    scriptFiles[i].close()

        