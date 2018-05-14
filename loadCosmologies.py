import numpy as np
import sys
from scipy import integrate
from scipy.interpolate import interp1d
from scipy.ndimage import gaussian_filter1d

def growthFunc(zz,HH):
    # this is the factor to convert the scale-less analytical growth factor to real
    H_zMax_scaleLess = (2.0/3.0*(1.0+zz[-1])**(3.0/2.0))
    ff = 1.0/(HH[-1] / H_zMax_scaleLess )**3 #This is essentially t_0^3 in the EdS formula, fitted at the highest z

    arg_EdS = lambda x: 1.0/(1.0+x)**(7.0/2.0)
    D_zMax_inf_EdS_raw, err = integrate.quad(arg_EdS, zz[-1], np.inf)
    D_zMax_inf_EdS_raw = 27.0/8.0 * D_zMax_inf_EdS_raw

    #print(err)
    
    D_zMax_inf_EdS = ff*D_zMax_inf_EdS_raw

    D_0_z = integrate.cumtrapz(np.divide(1.0+zz,HH**3),zz,initial=0)
    D_0_9 = integrate.trapz(np.divide(1.0+zz,HH**3),zz)

    D1 = HH/HH[0] * \
        (D_0_9 + D_zMax_inf_EdS - D_0_z) \
        / (D_0_9 + D_zMax_inf_EdS)
        
    return D1

times=[]
scaleFactors=[]
hubbleFactors=[]
omega_m_effs=[]

snapshotNames=[]

for sim in sims:
    globalSimParams=np.genfromtxt('../timeCurve_'+sim+'.txt',dtype='str')
    
    times.append(np.array([float(f) for f in globalSimParams[:,0]]))
    scaleFactors.append(np.array([float(f) for f in globalSimParams[:,1]]))
    hubbleFactors.append(np.array([float(f) for f in globalSimParams[:,3]]))
    omega_m_effs.append(np.array([float(f) for f in globalSimParams[:,4]]))

    snapshotNames.append(globalSimParams[:,2]) 
 
aVect=np.arange(0.104,1.0005,0.001)

for i in range(len(sims)):
    times[i]=interp1d(scaleFactors[i], times[i], kind='cubic', fill_value='extrapolate')(aVect)
    hubbleFactors[i]=interp1d(scaleFactors[i], hubbleFactors[i], kind='cubic', fill_value='extrapolate')(aVect)
    omega_m_effs[i]=interp1d(scaleFactors[i], omega_m_effs[i], kind='cubic', fill_value='extrapolate')(aVect)

    scaleFactors[i]=aVect

    
taus=[]
    
for i in range(len(sims)):
    
    taus.append(3*times[i][0]/scaleFactors[i][0]+integrate.cumtrapz(1.0/scaleFactors[i],times[i],initial=0))

D1_sim=[]
onePlusZD1_dz_sim=[]
onePlusZD1_dz_sim_smooth=[]

onePlusZD1_dtau_sim=[]
onePlusZD1_dtau_sim_smooth=[]

beta_sim=[]
beta_sim_smooth=[]
beta_sim_interp=[]

for i in range(len(sims)):
    
    D1_sim.append(np.flipud(growthFunc(np.flipud(1.0/scaleFactors[i]-1.0),np.flipud(hubbleFactors[i]))))
    
    onePlusZD1_dz_sim.append(np.gradient((1.0/scaleFactors[i])*D1_sim[i],edge_order=2)/
                             np.gradient(1.0/scaleFactors[i]-1.0,edge_order=2))
    
    onePlusZD1_dz_sim_smooth.append(gaussian_filter1d(onePlusZD1_dz_sim[i],sigma=10.0, mode='nearest'))
  
    onePlusZD1_dtau_sim.append(np.gradient((1.0/scaleFactors[i])*D1_sim[i],edge_order=2)/
                               np.gradient(taus[i],edge_order=2))
    
    onePlusZD1_dtau_sim_smooth.append(gaussian_filter1d(onePlusZD1_dtau_sim[i],sigma=10.0, mode='nearest'))

    beta_sim.append(np.gradient(np.log(D1_sim[i]),edge_order=2)/np.gradient(np.log(scaleFactors[i]),edge_order=2))
    
    beta_sim_smooth.append(gaussian_filter1d(beta_sim[i],sigma=10.0, mode='nearest'))
    
    beta_sim_interp.append(interp1d(taus[i], beta_sim_smooth[i], kind='cubic', fill_value='extrapolate'))

