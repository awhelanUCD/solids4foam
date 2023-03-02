import sys
import math
import os
import numpy as np
import copy

#height=30
#innerRadius=3
#notchedRadius=4
#outerRadius=5
#angle=5.0

height=20
innerRadius=5
notchedRadius=4
outerRadius=9
angle=2.5

combRadius=innerRadius+notchedRadius
innerValue=innerRadius/2.0
outerAngle=math.acos((innerRadius+notchedRadius-outerRadius)/notchedRadius)
print(math.degrees(outerAngle))

def rad(x):
    return x*(3.14/180)

mainStr='('+str(0)+' '+str(0)+' '+str(0)+')'+'\n'+'\n'

point_1='('+str((innerValue)*math.cos(rad(angle)))+' '+str(0)+' '+str(-(innerValue)*math.sin(rad(angle)))+')'+'\n'
point_2='('+str(innerRadius*math.cos(rad(angle)))+' '+str(0)+' '+str(-(innerRadius)*math.sin(rad(angle)))+')'+'\n'
point_3='('+str((combRadius-notchedRadius*math.cos(outerAngle/2.0))*math.cos(rad(angle)))+' '+str(notchedRadius*math.sin(outerAngle/2.0))+' '+str(-(combRadius-notchedRadius*math.cos(outerAngle/2.0))*math.sin(rad(angle)))+')'+'\n'
point_4='('+str((combRadius-(notchedRadius+innerValue)*math.cos(outerAngle/2.0))*math.cos(rad(angle)))+' '+str((notchedRadius+innerValue)*math.sin(outerAngle/2.0))+' '+str(-(combRadius-(notchedRadius+innerValue)*math.cos(outerAngle/2.0))*math.sin(rad(angle)))+')'+'\n'
point_5='('+str(outerRadius*math.cos(rad(angle)))+' '+str((notchedRadius*math.sin(outerAngle)))+' '+str(-outerRadius*math.sin(rad(angle)))+')'+'\n'
point_6='('+str(outerRadius*math.cos(rad(angle)))+' '+str(notchedRadius+innerValue)+' '+str(-outerRadius*math.sin(rad(angle)))+')'+'\n'
point_7='('+str(outerRadius*math.cos(rad(angle)))+' '+str(height)+' '+str(-outerRadius*math.sin(rad(angle)))+')'+'\n'
point_8='('+str((combRadius-(notchedRadius+innerValue)*math.cos(outerAngle/2.0))*math.cos(rad(angle)))+' '+str(height)+' '+str(-(combRadius-(notchedRadius+innerValue)*math.cos(outerAngle/2.0))*math.sin(rad(angle)))+')'+'\n'+'\n'

point_9='('+str(0)+' '+str(height)+' '+str(0)+')'+'\n'
point_10='('+str(0)+' '+str((notchedRadius+innerValue)*math.sin(outerAngle/2.0))+' '+str(0)+')'+'\n'+'\n'

point_11='('+str((innerValue)*math.cos(rad(angle)))+' '+str(0)+' '+str((innerValue)*math.sin(rad(angle)))+')'+'\n'
point_12='('+str(innerRadius*math.cos(rad(angle)))+' '+str(0)+' '+str((innerRadius)*math.sin(rad(angle)))+')'+'\n'
point_13='('+str((combRadius-notchedRadius*math.cos(outerAngle/2.0))*math.cos(rad(angle)))+' '+str(notchedRadius*math.sin(outerAngle/2.0))+' '+str((combRadius-notchedRadius*math.cos(outerAngle/2.0))*math.sin(rad(angle)))+')'+'\n'
point_14='('+str((combRadius-(notchedRadius+innerValue)*math.cos(outerAngle/2.0))*math.cos(rad(angle)))+' '+str((notchedRadius+innerValue)*math.sin(outerAngle/2.0))+' '+str((combRadius-(notchedRadius+innerValue)*math.cos(outerAngle/2.0))*math.sin(rad(angle)))+')'+'\n'
point_15='('+str(outerRadius*math.cos(rad(angle)))+' '+str((notchedRadius*math.sin(outerAngle)))+' '+str(outerRadius*math.sin(rad(angle)))+')'+'\n'
point_16='('+str(outerRadius*math.cos(rad(angle)))+' '+str(notchedRadius+innerValue)+' '+str(outerRadius*math.sin(rad(angle)))+')'+'\n'
point_17='('+str(outerRadius*math.cos(rad(angle)))+' '+str(height)+' '+str(outerRadius*math.sin(rad(angle)))+')'+'\n'
point_18='('+str((combRadius-(notchedRadius+innerValue)*math.cos(outerAngle/2.0))*math.cos(rad(angle)))+' '+str(height)+' '+str((combRadius-(notchedRadius+innerValue)*math.cos(outerAngle/2.0))*math.sin(rad(angle)))+')'+'\n'+'\n'

mainStr=mainStr+point_1+point_2+point_3+point_4+point_5+point_6+point_7+point_8+point_9+point_10+point_11+point_12+point_13+point_14+point_15+point_16+point_17+point_18

mainStr=mainStr+'Edges: ' +'\n'
arc1='arc 2 3 ('+str((combRadius-notchedRadius*math.cos(outerAngle/4.0))*math.cos(rad(angle)))+' '+str(notchedRadius*math.sin(outerAngle/4.0))+' '+str(-(combRadius-notchedRadius*math.cos(outerAngle/4.0))*math.sin(rad(angle)))+')'+'\n'
arc2='arc 12 13 ('+str((combRadius-notchedRadius*math.cos(outerAngle/4.0))*math.cos(rad(angle)))+' '+str(notchedRadius*math.sin(outerAngle/4.0))+' '+str((combRadius-notchedRadius*math.cos(outerAngle/4.0))*math.sin(rad(angle)))+')'+'\n'

arc3='arc 1 4 ('+str((combRadius-(notchedRadius+innerValue)*math.cos(outerAngle/3.0))*math.cos(rad(angle)))+' '+str((notchedRadius+innerValue)*math.sin(outerAngle/3.0))+' '+str(-(combRadius-(notchedRadius+innerValue)*math.cos(outerAngle/3.0))*math.sin(rad(angle)))+')'+'\n'
arc4='arc 11 14('+str((combRadius-(notchedRadius+innerValue)*math.cos(outerAngle/3.0))*math.cos(rad(angle)))+' '+str((notchedRadius+innerValue)*math.sin(outerAngle/3.0))+' '+str((combRadius-(notchedRadius+innerValue)*math.cos(outerAngle/3.0))*math.sin(rad(angle)))+')'+'\n'
arc5='arc 3 5 ('+str((combRadius-notchedRadius*math.cos(outerAngle/1.5))*math.cos(rad(angle)))+' '+str(notchedRadius*math.sin(outerAngle/1.5))+' '+str(-(combRadius-notchedRadius*math.cos(outerAngle/1.5))*math.sin(rad(angle)))+')'+'\n'
arc6='arc 13 15 ('+str((combRadius-notchedRadius*math.cos(outerAngle/1.5))*math.cos(rad(angle)))+' '+str(notchedRadius*math.sin(outerAngle/1.5))+' '+str((combRadius-notchedRadius*math.cos(outerAngle/1.5))*math.sin(rad(angle)))+')'+'\n'

arc7='arc 4 6 ('+str((combRadius-(notchedRadius+innerValue)*math.cos(outerAngle/1.5))*math.cos(rad(angle)))+' '+str((notchedRadius+innerValue)*math.sin(outerAngle/1.5))+' '+str(-(combRadius-(notchedRadius+innerValue)*math.cos(outerAngle/1.5))*math.sin(rad(angle)))+')'+'\n'
arc8='arc 14 16 ('+str((combRadius-(notchedRadius+innerValue)*math.cos(outerAngle/1.5))*math.cos(rad(angle)))+' '+str((notchedRadius+innerValue)*math.sin(outerAngle/1.5))+' '+str((combRadius-(notchedRadius+innerValue)*math.cos(outerAngle/1.5))*math.sin(rad(angle)))+')'+'\n'

mainStr=mainStr+arc1+arc2+arc3+arc4+arc5+arc6+arc7+arc8
print(mainStr)
text_file = open("partialBlockMesh", "w")
text_file.write(mainStr)
text_file.close()

