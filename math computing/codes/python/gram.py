import sys                                  
sys.path.insert(0, './CoordGeo')        
import numpy as np
import numpy.linalg as LA
import scipy.linalg as SA
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import pandas as pd
#local imports
from line.funcs import *
from triangle.funcs import *
from conics.funcs import circ_gen
#if using termux
import subprocess
import shlex
#end if
#Input parameters from excel file
df= pd.read_excel('vertices.xlsx')
#print(df ,"\n")
#Triangle Vertices
G_v= df.to_numpy()[:,:]
print("G_v=",G_v,"\n")
v1 = G_v[:, 0]
v2 = G_v[:, 1]
v3 = G_v[:, 2]
vn1=LA.norm(v1)
A=v1/vn1
B1= v2 - np.dot(v2.T, A) * A
vn2=LA.norm(B1)
B= B1/vn2
C1=v3-np.dot(v3.T , A)*A-np.dot(v3.T, B)*B
vn3=LA.norm(C1)
C= C1/vn3
#print(A)
#print(B)
#print(C)
print("***************gram schmidt matrix****************")
x=np.matrix((A,B,C))
print(x,"\n")
print("*********matrix after qr decomposition************")
y=x.T * x
print(y,"\n")
u1=round((LA.norm(y[:,0])))
u2=round((LA.norm(y[:,1])))
u3=round((LA.norm(y[:,2])))
#print(u1)
#print(u2)
#print(u3)
if(u1==u2==u3):
    print("A+B+C equally inclined to A,B,C")
else:
    print("not inclined equally")
