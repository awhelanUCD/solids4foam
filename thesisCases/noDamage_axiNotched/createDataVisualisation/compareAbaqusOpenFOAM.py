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
    plt.plot((time*1),(normalForce*72)/1000.0)



fileNames=['solidForcesup.dat']
plotData(fileNames)

text1=np.loadtxt('rf2abaqus.txt')
disp=text1[:,0]
forceR=text1[:,1]
plt.plot(disp,forceR/1e9,'--')

plt.ylabel('Force (kN)')
plt.xlabel('Displacement (mm)')
#plt.legend(fileNames)
plt.legend(['this work', 'Abaqus'])
plt.savefig('forceDisp.png',dpi=500)
plt.close()





