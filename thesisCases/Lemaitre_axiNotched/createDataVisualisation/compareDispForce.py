import math
import numpy as np
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
    plt.plot((time[:800]*1),(normalForce[:800])*72/1000)
    



fileNames1=['solidForcesup.dat']
plotData(fileNames1)

text1=np.loadtxt('abaqusDispForce.txt')
disp=text1[:,0]
force=text1[:,1]
plt.plot(disp[:800]*1,force[:800]/1e9,'--')

plt.xlabel('Displacement (mm)')
plt.ylabel('Force (kN)')
plt.legend(['OpenFOAM', 'Abaqus'])
plt.grid(linestyle='--',linewidth=0.5)
#plt.ylim([0, 3e8])

plt.savefig('1.1LemaitreCompare.png',dpi=500)








