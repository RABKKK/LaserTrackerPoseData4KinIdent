#Determine the relation between tracker and robot cooridnate frames
#Input: 
# WorldDef.xyz: measured pose as end-effector moves along coordinates of robot coordinate frames, used to compute orientation (pose data was grabbed from laser tracker using Metrolog software and saved in the .xyz format)
# RobPoseatRef.csv: Pose of robot flange pose (position and orientation) from robot coontroller
# TransToolPosit.csv: The position of flange estimated during the identification of relation between tool and flange frames
#Output
#TransTractoWorld.csv: Transformation from tracker to robot coordinate frame
import re
import csv
from numpy import savetxt
import cv2 as cv
import numpy as np
def rotmat(th,cond):
    if cond=='x':
        thr=th*np.pi/180
        R=np.array([[1,0,0],[0,np.math.cos(thr),-np.math.sin(thr)],[0,np.math.sin(thr),np.math.cos(thr)]])
    elif cond=='y':
        thr=th*np.pi/180
        R=np.array([[np.math.cos(thr),0,np.math.sin(thr)],[0,1,0],[-np.math.sin(thr),0,np.math.cos(thr)]])
    elif cond=='z':
        thr=th*np.pi/180
        R=np.array([[np.math.cos(thr),-np.math.sin(thr),0],[np.math.sin(thr),np.math.cos(thr),0],[0,0,1],])
    else:
        print('Input error: enter x for x axis rot, y for y axis etc')
        R=np.eye(3,3)
    return(R)

def transmat(rotpar):
    print('Fanuc x, y, z')
    R1=rotmat(rotpar[3],'x')
    R2=rotmat(rotpar[4],'y')
    R3=rotmat(rotpar[5],'z')
    R=np.matmul(R3,np.matmul(R2,R1))
    T=np.concatenate((R,np.array([[rotpar[0]],[rotpar[1]],[rotpar[2]]])),1)
    Th=np.vstack((T,np.array([0,0,0,1])))
    return(Th)
def nestposetotransfmat(N1,N2,N3):
    #N2 origin, N2N1 x axis, N2N3 on xy plane
    xax=(N1-N2)/np.linalg.norm(N1-N2)
    xyplvec=(N3-N2)/np.linalg.norm(N3-N2)
    zax=np.cross(xax,xyplvec)/np.linalg.norm(np.cross(xax,xyplvec))
    yax=np.cross(zax,xax)/np.linalg.norm(np.cross(zax,xax))
    Thmat=np.hstack((np.transpose(np.vstack((np.vstack((xax,yax)),zax))),np.array([[N2[0]],[N2[1]],[N2[2]]])))
    Thmat=np.vstack((Thmat,np.array([0,0,0,1])))
    print('Transform mat',Thmat)
    return(Thmat)
def toolposit(T1,T2,T3,T4):
    matpx=np.vstack((np.vstack((np.hstack((T1[0,:3],np.array([-1,0,0]))),np.hstack((T2[0,:3],np.array([-1,0,0]))))),np.hstack((T3[0,:3],np.array([-1,0,0])))))
    matpx=np.vstack((matpx,np.hstack((T4[0,:3],np.array([-1,0,0])))))
    matpy=np.vstack((np.vstack((np.hstack((T1[1,:3],np.array([0,-1,0]))),np.hstack((T2[1,:3],np.array([0,-1,0]))))),np.hstack((T3[1,:3],np.array([0,-1,0])))))
    matpy=np.vstack((matpy,np.hstack((T4[1,:3],np.array([0,-1,0])))))
    matpz=np.vstack((np.vstack((np.hstack((T1[2,:3],np.array([0,0,-1]))),np.hstack((T2[2,:3],np.array([0,0,-1]))))),np.hstack((T3[2,:3],np.array([0,0,-1])))))
    matpz=np.vstack((matpz,np.hstack((T4[2,:3],np.array([0,0,-1])))))	
    cpx=-np.array([[T1[0,3]],[T2[0,3]],[T3[0,3]],[T4[0,3]]])
    cpy=-np.array([[T1[1,3]],[T2[1,3]],[T3[1,3]],[T4[1,3]]])
    cpz=-np.array([[T1[2,3]],[T2[2,3]],[T3[2,3]],[T4[2,3]]])
    cmatpxyz=np.vstack((matpx,matpy))
    cmatpxyz=np.vstack((cmatpxyz,matpz))
    ccpxyz=np.vstack((cpx,cpy))
    ccpxyz=np.vstack((ccpxyz,cpz))
    txyzpxyz=np.matmul(np.linalg.pinv(cmatpxyz),ccpxyz)
    return(txyzpxyz[:3],txyzpxyz[3:])

def calctooldeferr(ThMatP1P4,transvec,toolposit):
    cnt=0
    terr=np.zeros(int(ThMatP1P4.shape[0]/4))
    ThNesttoTool=np.eye(4,4)
    ThNesttoTool[0,3]=transvec[0]
    ThNesttoTool[1,3]=transvec[1]
    ThNesttoTool[2,3]=transvec[2]
    for line in ThMatP1P4:
        cnt+=1
        if cnt%4==0:
              t1=np.matmul(ThMatP1P4[cnt-4:cnt,:],ThNesttoTool)
              print(t1[:3,3],toolposit[:,0])
              terr[int(cnt/4)-1]=np.linalg.norm(t1[:3,3]-toolposit[:,0])
    return(terr) 

#Orientation def
fnamorient='WorldDef.xyz'
OrientPosVal=np.zeros((4,3))
with open(fnamorient) as f:
    cnt=0
    for line in f:
        cnt+=1
        spltcont=re.split('[;]|[;;]|[;;;]|[;;;;]|[;;;;;]|[\n]',line)
        OrientPosVal[cnt-1,:3]=np.array([float(spltcont[1]),float(spltcont[2]), float(spltcont[3])])

xax=(OrientPosVal[1,:]-OrientPosVal[0,:])/np.linalg.norm(OrientPosVal[1,:]-OrientPosVal[0,:])#Obtained by moving along x axis of tool, value from metrolog, row 0 to row 1 is x axis
print('X axis vector',xax)
zax=(OrientPosVal[3,:]-OrientPosVal[2,:])/np.linalg.norm(OrientPosVal[3,:]-OrientPosVal[2,:])##Obtained by moving along z axis of tool, value from metrolog, row 0 to row 2 is z axis
print('z axis vector',zax)
yax=np.cross(zax,xax)/np.linalg.norm(np.cross(zax,xax))
zax=np.cross(xax,yax)/np.linalg.norm(np.cross(xax,yax))
RmatTractoWorld=np.transpose(np.vstack((np.vstack((xax,yax)),zax)))
#Considering pose 4
robpose=np.loadtxt('RobPoseatRef.csv')#robot pose at [4] according to controller
toolPosit=np.loadtxt('TransToolPosit.csv')
ThWorldtoFlange=transmat(robpose)
transTractoWorld=np.matmul(RmatTractoWorld,-ThWorldtoFlange[:3,3])+toolPosit[1,:]
print('Tracker to World rotation',RmatTractoWorld,'Translation',transTractoWorld)
ThTractoWorld=np.hstack((RmatTractoWorld,np.array([[transTractoWorld[0]],[transTractoWorld[1]],[transTractoWorld[2]]])))
print('Tractor to World coord',ThTractoWorld)
np.savetxt('TransTractoWorld.csv',ThTractoWorld)

    

