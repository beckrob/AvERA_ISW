import numpy as np
import random
import datetime
import math
import sys
from scipy import integrate
from scipy.interpolate import interp1d
from joblib import Parallel, delayed
import healpy as hp


sims=[sys.argv[1]]
startID=int(sys.argv[2])

simulationType=int(sys.argv[3])

outRedshiftList=[0.5,1.0,1.5,2.0,2.5,3.0,4.0,5.0,6.0,7.0,8.0]

exec(open('./loadCosmologies.py').read())
    
milleniumFiles=np.loadtxt('../MilleniumXXL/zlist_da.txt')
   
startIndex=14
endIndex=64
   
if simulationType==0:

    inPath= './LC_results/'
    outPath= './ISW_map/'
    
elif simulationType==1 or simulationType==2:
    
    if simulationType==1:
    
        inPath= './LC_results_XXL/'
        outPath= './ISW_map/XXL_'
        
    elif simulationType==2:
    
        inPath= './LC_results_XXL_lin/'
        outPath= './ISW_map_lin/XXL_'
    
    #tau_a_func=interp1d(scaleFactors[0], taus[0], kind='linear', fill_value='extrapolate')
    #H_tau_func=interp1d(taus[0], hubbleFactors[0], kind='linear', fill_value='extrapolate')
    
    #h=H_tau_func(tau_a_func(milleniumFiles[endIndex-1,1]))/100.0
    
    h_MXXL=0.73
    
    boxSize=3000/h_MXXL
    
    M_pointmass=8.456e9 #*0.73/h #In solar masses
    
    pointMassCount=303464448000.0
    
    
a_tau_func=interp1d(taus[0], scaleFactors[0], kind='linear', fill_value='extrapolate')
    
healPixResolution=64
startPointNum = 3

(thetas,phis)=hp.pix2ang(healPixResolution, np.arange(hp.nside2npix(healPixResolution)))

#Physical units in SI
G_const=6.6740831e-11
M_Sun=1.98892e30


#Mass_LCDM=2.024938
#Mass BR=2.029469

#Parameter choices adopted from cosmolopy (https://roban.github.io/CosmoloPy/docAPI/cosmolopy.constants-module.html)
speedOfLightSI = 299792458.0 
MpcInMeters = 3.08568025e22
GyrInSeconds = 3.1536e16

speedOfLightMpcGyr = speedOfLightSI/MpcInMeters*GyrInSeconds

resultSImultiplier=1000*MpcInMeters

CMBTemp=2.726
        
for i in range(len(sims)):
    
    lightCurve=None
    massPrefactor=0.0
    
    if simulationType==0:
    
        for j in range(len(snapshotNames[i])):
        
            lightCurvePart=np.load(inPath+sims[i]+'/'+snapshotNames[i][j][:-4]+'_'+str(startID)+'.npy')
            
            if lightCurve is None:
                
                lightCurve=lightCurvePart
                
            else:
                
                lightCurve=np.vstack((lightCurve,lightCurvePart))

        H0_SI=hubbleFactors[i][-1]*1000/MpcInMeters
        
        omega_m=omega_m_effs[i][-1]                
                
        massPrefactor=3.0*(H0_SI**2)*omega_m
        
    elif simulationType==1 or simulationType==2:
        
        for j in range(startIndex,endIndex):
            
            lightCurvePart=np.load(inPath+sims[i]+'/density_1024_0'+str(j)+'_'+str(startID)+'.npy')
            
            if lightCurve is None:
                
                lightCurve=lightCurvePart
                
            else:
                
                lightCurve=np.vstack((lightCurve,lightCurvePart))
        
        massPrefactor=8.0*math.pi*G_const*pointMassCount*M_pointmass*M_Sun/(boxSize*MpcInMeters*milleniumFiles[endIndex-1,1])**3

    if simulationType==0:
    
        deltaT_matrix=[]
        
    deltaT_linear_matrix=[]
    
    for dirID in range(len(thetas)):
    
        theta=thetas[dirID]
        phi=phis[dirID]

        singleLightCurve=np.flipud(lightCurve[np.logical_and(lightCurve[:,0]==theta,lightCurve[:,1]==phi),:])

        dtau_SI=singleLightCurve[:,2]*GyrInSeconds
        
        redshiftCurve=1.0/a_tau_func(singleLightCurve[:,2])-1.0
        
        it=0
        outIndices=[]

        for zId in range(len(redshiftCurve)):
             
            if it<len(outRedshiftList) and redshiftCurve[zId]>outRedshiftList[it]:
                
                outIndices.append(zId)
                it+=1
                
        outIndices.append(len(redshiftCurve)-1)    
        
        if simulationType==0:
            
            dphi_dtau_SI=singleLightCurve[:,6]*resultSImultiplier*(-1.0)/(speedOfLightSI**2)*massPrefactor
            #This includes all constant multipliers
        
            deltaT=integrate.cumtrapz(dphi_dtau_SI,dtau_SI,initial=0)*CMBTemp
            
            deltaT_matrix.append([theta,phi]+list(deltaT[outIndices]))
        
        
        dphi_dtau_SI_linear=singleLightCurve[:,7]*resultSImultiplier*(-1.0)/(speedOfLightSI**2)*massPrefactor
    
        #WARNING: now the curves are expected to be multiplied by 1-beta beforehand
        #dphi_dtau_SI_linear*=(1.0-beta_sim_interp[i](singleLightCurve[:,2]))
    
        deltaT_linear=integrate.cumtrapz(dphi_dtau_SI_linear,dtau_SI,initial=0)*CMBTemp
        
        deltaT_linear_matrix.append([theta,phi]+list(deltaT_linear[outIndices]))
        
    if simulationType==0:
        deltaT_matrix=np.array(deltaT_matrix)
        np.savetxt(outPath+'deltaT_'+sims[i]+'_'+str(startID)+'.txt',deltaT_matrix)
                
    
    deltaT_linear_matrix=np.array(deltaT_linear_matrix)
    np.savetxt(outPath+'deltaT_linear_'+sims[i]+'_'+str(startID)+'.txt',deltaT_linear_matrix)

