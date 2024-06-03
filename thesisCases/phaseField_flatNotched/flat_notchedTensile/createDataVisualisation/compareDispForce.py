##import sys
##import os
#import numpy as np
##import copy
##import random



##lhs(n, [samples, criterion, iterations])
#x=lhs(2 , samples=5,criterion="m",iterations=40)
#print x

import sys
import math
import os
import numpy as np
import copy
import pandas as pd
import matplotlib.pyplot as plt

def plotData(fileNames):
  area=20.4768
  for file in fileNames:
    forceData=np.loadtxt(file,usecols=(0,1,2,3,4))
    print(forceData[0])
    xForce=[]
    normalForce=[]
    time=[]
    yForce=[]
    for data in forceData:
        time.append(data[0])
        xForce.append(data[1])
        normalForce.append(data[4])
        yForce.append(data[2])

    xForce=np.array(xForce)
    yForce=np.array(yForce)
    normalForce=np.array(normalForce)
    time=np.array(time)
    totalForce=(normalForce**2+xForce**2)**0.5
    print(totalForce)
   
    plt.plot(time[:306]*1.143/76.2,normalForce[:306]/(area))
  
    


#text1=np.loadtxt('0.01abaqus.txt')
#disp=text1[:,0]
#force=text1[:,1]
#plt.plot(disp,force/1e9,'--')


fileNames1=['solidForcesup.dat']#,'_solidForcesup.dat','sigma_solidForcesup.dat']
plotData(fileNames1)

df=pd.read_csv('bordenData.csv')
strain=df['Strain']
stress=df['Stress']
plt.plot(strain,stress,'--')#'o',ms=2)

df=pd.read_csv('eldahshanData.csv')
strain=df['Strain']
stress=df['Stress']
plt.plot(strain,stress,'o',ms=2)

plt.xlabel('Normalised Strain')
plt.ylabel('Normalised Stress')
plt.legend(['OpenFOAM','Borden et Al. (2016)','Eldahshan et Al. (2021)'])
plt.grid(linestyle='--',linewidth=0.5)
#plt.ylim([0, 3e8])

plt.savefig('phaseCompare.png',dpi=200)








