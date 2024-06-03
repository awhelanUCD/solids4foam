import sys
import math
import os
import numpy as np
import copy

height=76.2
innerRadius=8.64
notchedRadius=4.06
outerRadius=12.7
angle=5.07
thickness=2.0

combRadius=innerRadius+notchedRadius
innerValue=innerRadius/2.0
outerAngle=math.acos((innerRadius+notchedRadius-outerRadius)/notchedRadius)
print(math.degrees(outerAngle))

def rad(x):
    return x*(3.14/180)

mainStr='('+str(0)+' '+str(0)+' '+str(0)+')'+'\n'+'\n'

point_1='('+str(innerValue)+' '+str(0)+' '+str(0)+')'+'\n'
point_2='('+str(innerRadius)+' '+str(0)+' '+str(0)+')'+'\n'
point_3='('+str((combRadius-notchedRadius*math.cos(outerAngle/2.0)))+' '+str(notchedRadius*math.sin(outerAngle/2.0))+' '+str(0)+')'+'\n'
point_4='('+str((combRadius-(notchedRadius+innerValue)*math.cos(outerAngle/2.0)))+' '+str((notchedRadius+innerValue)*math.sin(outerAngle/2.0))+' '+str(0)+')'+'\n'
point_5='('+str(outerRadius)+' '+str((notchedRadius*math.sin(outerAngle)))+' '+str(0)+')'+'\n'
point_6='('+str(outerRadius)+' '+str(notchedRadius+innerValue)+' '+str(0)+')'+'\n'
point_7='('+str(outerRadius*math.cos(rad(angle)))+' '+str(height)+' '+str(0)+')'+'\n'
point_8='('+str((combRadius-(notchedRadius+innerValue)*math.cos(outerAngle/2.0)))+' '+str(height)+' '+str(0)+')'+'\n'+'\n'

point_9='('+str(0)+' '+str(height)+' '+str(0)+')'+'\n'
point_10='('+str(0)+' '+str((notchedRadius+innerValue)*math.sin(outerAngle/2.0))+' '+str(0)+')'+'\n'+'\n'

point_11='('+str(0)+' '+str(0)+' '+str(thickness)+')'+'\n'+'\n'

point_12='('+str(innerValue)+' '+str(0)+' '+str(thickness)+')'+'\n'
point_13='('+str(innerRadius)+' '+str(0)+' '+str(thickness)+')'+'\n'
point_14='('+str((combRadius-notchedRadius*math.cos(outerAngle/2.0)))+' '+str(notchedRadius*math.sin(outerAngle/2.0))+' '+str(thickness)+')'+'\n'
point_15='('+str((combRadius-(notchedRadius+innerValue)*math.cos(outerAngle/2.0)))+' '+str((notchedRadius+innerValue)*math.sin(outerAngle/2.0))+' '+str(thickness)+')'+'\n'
point_16='('+str(outerRadius)+' '+str((notchedRadius*math.sin(outerAngle)))+' '+str(thickness)+')'+'\n'

point_17='('+str(outerRadius)+' '+str(notchedRadius+innerValue)+' '+str(thickness)+')'+'\n'
point_18='('+str(outerRadius*math.cos(rad(angle)))+' '+str(height)+' '+str(thickness)+')'+'\n'
point_19='('+str((combRadius-(notchedRadius+innerValue)*math.cos(outerAngle/2.0)))+' '+str(height)+' '+str(thickness)+')'+'\n'+'\n'

point_20='('+str(0)+' '+str(height)+' '+str(thickness)+')'+'\n'
point_21='('+str(0)+' '+str((notchedRadius+innerValue)*math.sin(outerAngle/2.0))+' '+str(thickness)+')'+'\n'+'\n'

mainStr=mainStr+point_1+point_2+point_3+point_4+point_5+point_6+point_7+point_8 \
+point_9+point_10+point_11+point_12+point_13+point_14+point_15+point_16+point_17+point_18 \
+point_19+point_20+point_21

mainStr=mainStr+'Edges: ' +'\n'
arc1='arc 2 3 ('+str((combRadius-notchedRadius*math.cos(outerAngle/4.0)))+' '+str(notchedRadius*math.sin(outerAngle/4.0))+' '+str(0)+')'+'\n'
arc2='arc 13 14 ('+str((combRadius-notchedRadius*math.cos(outerAngle/4.0)))+' '+str(notchedRadius*math.sin(outerAngle/4.0))+' '+str(thickness)+')'+'\n'

arc3='arc 1 4 ('+str((combRadius-(notchedRadius+innerValue)*math.cos(outerAngle/3.0)))+' '+str((notchedRadius+innerValue)*math.sin(outerAngle/3.0))+' '+str(0)+')'+'\n'
arc4='arc 12 15 ('+str((combRadius-(notchedRadius+innerValue)*math.cos(outerAngle/3.0)))+' '+str((notchedRadius+innerValue)*math.sin(outerAngle/3.0))+' '+str(thickness)+')'+'\n'

arc5='arc 3 5 ('+str((combRadius-notchedRadius*math.cos(outerAngle/1.5)))+' '+str(notchedRadius*math.sin(outerAngle/1.5))+' '+str(0)+')'+'\n'
arc6='arc 14 16 ('+str((combRadius-notchedRadius*math.cos(outerAngle/1.5)))+' '+str(notchedRadius*math.sin(outerAngle/1.5))+' '+str(thickness)+')'+'\n'

arc7='arc 4 6 ('+str((combRadius-(notchedRadius+innerValue)*math.cos(outerAngle/1.5)))+' '+str((notchedRadius+innerValue)*math.sin(outerAngle/1.5))+' '+str(0)+')'+'\n'
arc8='arc 15 17 ('+str((combRadius-(notchedRadius+innerValue)*math.cos(outerAngle/1.5)))+' '+str((notchedRadius+innerValue)*math.sin(outerAngle/1.5))+' '+str(thickness)+')'+'\n'

mainStr=mainStr+arc1+arc2+arc3+arc4+arc5+arc6+arc7+arc8
print(mainStr)
text_file = open("partialBlockMesh", "w")
text_file.write(mainStr)
text_file.close()

