import numpy as np
import random
import math
import sys
import numpy.fft as fft
import healpy as hp
from scipy import integrate
from scipy.interpolate import interp1d

def loadMilleniumDensity(filePath):

    datafile=open(filePath, 'rb')

    scaleFactor=np.fromfile(datafile, dtype=np.dtype('f4'), count=1, sep='')
    cellNum=np.fromfile(datafile, dtype=np.dtype('i4'), count=1, sep='')

    density=8.0*np.fromfile(datafile, dtype=np.dtype('f4'), count=cellNum[0]**3, sep='').reshape(([cellNum[0]]*3))

    datafile.close()

    return density

def interpolateValueOnGrid3D(scaledCoord, nullIndex, gridValues, gridSize):

    value=0.0

    for i in range(2):
        for j in range(2):
            for k in range(2):
                
                value+=np.prod(np.abs(1-scaledCoord-np.array([i,j,k])))*gridValues[(nullIndex[0]+i)%gridSize,(nullIndex[1]+j)%gridSize,(nullIndex[2]+k)%gridSize]
                

    return value  
    
sims=[sys.argv[7]]
    
gridDensities=loadMilleniumDensity(sys.argv[1]) # In count/cell volume

exec(open('./loadCosmologies.py').read())


H_tau_func=interp1d(taus[0], hubbleFactors[0], kind='cubic', fill_value='extrapolate')
D_tau_func=interp1d(taus[0], D1_sim[0], kind='cubic', fill_value='extrapolate')


tauStart=float(sys.argv[2]) #In Gyr
tauEnd=float(sys.argv[3]) #In Gyr

tauInitial=float(sys.argv[4]) #In Gyr
tauFinal=float(sys.argv[5]) #In Gyr

h_MXXL=0.73

#boxSize=4109.6 # size of the box in Mpc
#boxSize=3000/(H_tau_func(tauFinal)/100.0) # size of the box in Mpc

boxSize=3000/h_MXXL # scaling box size in Mpc to 

pointMassCount=303464448000.0

healPixResolution=64
startPointNum = 3

maxStepSize=0.75 # in Mpc

dim=3

#Parameter choices adopted from cosmolopy (https://roban.github.io/CosmoloPy/docAPI/cosmolopy.constants-module.html)
speedOfLightSI = 299792458.0 
MpcInMeters = 3.08568025e22
GyrInSeconds = 3.1536e16

linearGridSize=gridDensities.shape[0]

speedOfLightMpcGyr = speedOfLightSI/MpcInMeters*GyrInSeconds

lightTravelDistance1=speedOfLightMpcGyr*(tauFinal-tauStart)
lightTravelDistance2=speedOfLightMpcGyr*(tauFinal-tauEnd)

random.seed(2+5*3*45+3452*11*17+31)
for k in range(3000):
    dummy=random.random()


#Generate random starting points  
startPoints = np.empty((startPointNum,dim))
for i in range(startPointNum):
    for j in range(dim):
        startPoints[i,j]=random.random()*boxSize

        
(thetas,phis)=hp.pix2ang(healPixResolution, np.arange(hp.nside2npix(healPixResolution)))

averageDensity=pointMassCount/(linearGridSize)**3 # In count/cell volume

delta_r=gridDensities/averageDensity-1


k_vect=2*math.pi*fft.fftfreq(linearGridSize,boxSize/linearGridSize) #This is in units of 1/Mpc
k_vectLast=2*math.pi*fft.rfftfreq(linearGridSize,boxSize/linearGridSize)


delta_k=fft.rfftn(delta_r)
   

k2_inverse=np.empty(delta_k.shape, dtype=np.complex_)

if (dim==3):

    for i in range(delta_k.shape[0]):
        for j in range(delta_k.shape[1]):
            for k in range(delta_k.shape[2]):
                
                k2_inverse[i,j,k]=1.0/(k_vect[i]*k_vect[i]+
                                       k_vect[j]*k_vect[j]+
                                       k_vectLast[k]*k_vectLast[k])
                
    k2_inverse[0,0,0]=0.0

if (dim==2):

    for i in range(delta_k.shape[0]):
        for j in range(delta_k.shape[1]):
            
            k2_inverse[i,j]=1.0/(k_vect[i]*k_vect[i]+
                                 k_vectLast[j]*k_vectLast[j])

    k2_inverse[0,0]=0.0
            
            

phiDot_k_linear = k2_inverse*delta_k

phiDot_r_linear = fft.irfftn(phiDot_k_linear)


for startID in range(startPointNum):

    sampleResults=[]
    
    for dirID in range(len(thetas)):
    
        theta=thetas[dirID]
        phi=phis[dirID]
        
        sampleFrom=np.empty(3)
        
        sampleFrom[0]=startPoints[startID,0]+math.sin(theta)*math.cos(phi)*lightTravelDistance1
        sampleFrom[1]=startPoints[startID,1]+math.sin(theta)*math.sin(phi)*lightTravelDistance1
        sampleFrom[2]=startPoints[startID,2]+math.cos(theta)*lightTravelDistance1


        sampleTo=np.empty(3)

        sampleTo[0]=startPoints[startID,0]+math.sin(theta)*math.cos(phi)*lightTravelDistance2
        sampleTo[1]=startPoints[startID,1]+math.sin(theta)*math.sin(phi)*lightTravelDistance2
        sampleTo[2]=startPoints[startID,2]+math.cos(theta)*lightTravelDistance2
        
        drTotal=(sampleTo-sampleFrom)
        distTotal=math.sqrt(np.sum(drTotal*drTotal))

        interpolationSteps=int(distTotal/maxStepSize)+1

        for s in range(interpolationSteps):
            
            tau=tauStart+(tauEnd-tauStart)*s/interpolationSteps
            
            sample=sampleFrom+drTotal*s/interpolationSteps
            
            sampleShifted=sample-np.array([math.floor(f) for f in sample/boxSize])*boxSize
            
            cubeIndices=np.array([int(math.floor(f)) for f in (sampleShifted/(boxSize/linearGridSize)-0.5)])
            
            scaledCoords=sampleShifted/(boxSize/linearGridSize)-(cubeIndices+0.5)
            
            #Scaling overdensity (and thus phiDot) by the growth factor, normalized to the time of the initial snapshot
            phiDot_result_linear=H_tau_func(tau)*interpolateValueOnGrid3D(scaledCoords,cubeIndices,phiDot_r_linear,linearGridSize)*(1.0-beta_sim_interp[0](tau))*D_tau_func(tau)/D_tau_func(tauInitial)
            
            sampleResults.append([theta,phi,tau,sample[0],sample[1],sample[2],-9999.0,phiDot_result_linear])
    

    np.save(sys.argv[6]+str(startID)+'.npy', np.array(sampleResults))

