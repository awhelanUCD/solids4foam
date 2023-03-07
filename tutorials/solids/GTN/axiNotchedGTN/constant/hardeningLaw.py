import sys
import math
import os
import numpy as np
import copy
import pandas as pd
import matplotlib.pyplot as plt

#global yinf=YINF
#global y0=Y0
#global beta=BETA

def hardeningLaw(x):
    y1=1009e6*(0.003552+x)**0.1485
    y2=515.5e6+274.9e6*(1-2.71828**(-13.75*x))
    y=0.7*y1+0.3*y2
    y=589e6*(0.0001+x)**0.216
   # return '('+str(x)+'    '+str(y)+')'
    return '('+str(x)+'    '+str(y)+')'

strainValues=[0.0001*(i) for i in range(30000)]

#print(strainValues)
mainStr='( \n'
for i in strainValues:
    str1=hardeningLaw(i)
    mainStr=mainStr+str1+'\n'

mainStr=mainStr+')'
#print(mainStr)
text_file = open("plasticStrainVsYieldStress", "w")
text_file.write(mainStr)
text_file.close()

