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
  for file in fileNames:
    forceData=np.loadtxt(file,usecols=(0,1,2,3,4))
    print(forceData[0])
    xForce=[]
    normalForce=[]
    time=[]
    for data in forceData:
        time.append(data[0])
        xForce.append(data[1])
        normalForce.append(data[4])

    xForce=np.array(xForce)
    normalForce=np.array(normalForce)
    time=np.array(time)
    totalForce=(normalForce**2+xForce**2)**0.5
    print(totalForce)
    plt.plot((time*1),(normalForce*72)/1000.0)
    
text1=np.loadtxt('abaqus.txt')
disp=text1[:,0]
force=text1[:,1]
plt.plot(disp,force/1e9,'--')




fileNames1=['solidForcesup.dat']#,'_solidForcesup.dat','sigma_solidForcesup.dat']
plotData(fileNames1)

plt.xlabel('Displacement (mm)')
plt.ylabel('Force (kN)')
plt.legend(['Abaqus','OpenFOAM'])
plt.grid(linestyle='--',linewidth=0.5)
plt.xlim([0, 0.5])

plt.savefig('forceDispGTN.png',dpi=200)








