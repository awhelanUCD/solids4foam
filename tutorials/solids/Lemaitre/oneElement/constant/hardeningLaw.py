import math

def hardeningLaw(x):

    y=200e6+10e9*x
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

