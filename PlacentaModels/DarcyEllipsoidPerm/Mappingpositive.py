import numpy as np
from parameter import *
import time
import sys
tic = time.time()


ConEl=np.loadtxt('ConductivityEachEl_2ndsmooth.txt')
PEl=np.loadtxt('porosityEachEl.txt')
P=np.loadtxt('nodeSpacialData.txt')
P=P[:,1:4]*10# convert cm to mm
Px=P[:,0]
Py=P[:,1]
Pz=P[:,2]

halfPx= (abs(np.amin(Px))+abs(np.amax(Px)))/2
halfPy= (abs(np.amin(Py))+abs(np.amax(Py)))/2
halfPz= (abs(np.amin(Pz))+abs(np.amax(Pz)))/2

Px=Px+halfPx 
Py=Py+halfPy
Pz=Pz+halfPz

newP=np.vstack((Px,Py,Pz)).T
f=open("mappedK&P_final.txt",'w')
errorcount=0
for j in range(0,len(newP)):


  Endpoints = newP[j]

  if Endpoints[0]>nel_x or Endpoints[1]>nel_y or Endpoints[2]>nel_z:
     errorcount=errorcount+1
  
  num_x = np.floor((Endpoints[0]-x_min)/x_width)+1;
  num_y = np.floor((Endpoints[1]-y_min)/y_width)+1;
  num_z = np.floor((Endpoints[2]-z_min)/z_width)+1;
   
  if num_x > nel_x:
     num_x = nel_x

  if num_y >= nel_y:
     num_y = nel_y
   
  if num_z >= nel_z:
     num_z = nel_z

  T_e = ((num_z-1)*nel_x*nel_y + (num_y-1)*nel_x + num_x)-1
  E_found=T_e+1
  f.write("%s %s %s %s %s %s \n" %(j+1,P[j,0]/10,P[j,1]/10,P[j,2]/10,ConEl[T_e],PEl[T_e]))
  
f.close()
print 'errorcount', errorcount

'''
KK=np.loadtxt('mappedK&P_final.txt')

f=open("MappedK.exnode",'w')

f.write (" Group name: mapped_Conductivity\n")
f.write (" #Fields=2\n")
f.write (" 1) coordinates, coordinate, rectangular cartesian, #Components=3\n")
f.write (" x.  Value index=1, #Derivatives=0\n")
f.write (" y.  Value index=2, #Derivatives=0\n")
f.write (" z.  Value index=3, #Derivatives=0\n")
f.write (" 2) general, field, rectangular cartesian, #Components=1\n")
f.write ("  1.  Value index=4, #Derivatives=0\n")

for ee in range(0,len(P)):
    f.write("Node:  "        "%s\n" %(ee+1))
    f.write("          %s\n" %P[ee,0])
    f.write("          %s\n" %P[ee,1])
    f.write("          %s\n" %P[ee,2])
    f.write("          %s\n" %KK[ee,4])




f.close()


f=open("MappedP.exnode",'w')# this is k value in each node of mesh

f.write (" Group name: mapped_Porosity\n")
f.write (" #Fields=2\n")
f.write (" 1) coordinates, coordinate, rectangular cartesian, #Components=3\n")
f.write (" x.  Value index=1, #Derivatives=0\n")
f.write (" y.  Value index=2, #Derivatives=0\n")
f.write (" z.  Value index=3, #Derivatives=0\n")
f.write (" 2) general, field, rectangular cartesian, #Components=1\n")
f.write ("  1.  Value index=4, #Derivatives=0\n")

for ee in range(0,len(P)):
    f.write("Node:  "        "%s\n" %(ee+1))
    f.write("          %s\n" %P[ee,0])
    f.write("          %s\n" %P[ee,1])
    f.write("          %s\n" %P[ee,2])
    f.write("          %s\n" %KK[ee,5])




f.close()
'''


toc = time.time() - tic
print 'elapsed_time',toc



