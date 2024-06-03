import math
import numpy as np
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
    print((totalForce*72)/1000.0)
    plt.plot((time*1.143),(normalForce)/1000.0)



fileNames=['solidForcesup.dat']#,'bathe_solidForcesup.dat']#,'sol_solidForcesup.dat']

plotData(fileNames)


text1=np.loadtxt('rf2abaqus.txt')
disp=text1[:,0]
forceR=text1[:,1]
plt.plot(disp*1.143,forceR/1e9,'--')

plt.ylabel('Force (kN)')
plt.xlabel('Displacement (mm)')
#plt.legend(fileNames)
plt.legend(['this work','Abaqus'])
#plt.legend(['OpenFOAM','OpenFOAMNeoHookean','kirchoff','bathe Bekaert-v2.0','Abaqus'])
plt.savefig('forceDisp.png',dpi=500)
plt.close()





