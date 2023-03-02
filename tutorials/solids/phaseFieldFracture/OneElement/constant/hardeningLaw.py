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
    y=320e6+688e6*x
    #y=345+7148*x
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

