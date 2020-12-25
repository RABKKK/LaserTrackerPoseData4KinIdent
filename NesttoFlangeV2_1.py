#estimation of relation between Nest to Flange coordinate frames
#Input:
#ToolOrientDef.xyz: measured after by moving the end-effector fixing the position of flange origin and ToolTransPose.xyz was measured after moving the flange along its own x and z coordinate axes. (pose data was grabbed from laser tracker using Metrolog software and saved in the .xyz format)
#Output TransFlangtoNest.csv and TransToolPosit.csv (used for estimating transformation from robot to tracker coordinate frames)
#Use in conjunction with article R. A. Boby, A. Klimchik, Combination of Geometric and Parametric Approaches for Kinematic Identification of an Industrial Robot, RCIM
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
    #print('Transform mat',Thmat)
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

#Main function

NoNest=4#Number of nests positions measured
fnam='ToolTransPose.xyz'#Measured positions of Nest
NestVal=np.zeros((16,3))#Assuming 4 nests on end-effector, to cater to occlusion
ThMatP1P4=np.zeros((16,4))
with open(fnam) as f:
    cnt=0
    for line in f:
        cnt+=1
        spltcont=re.split('[;]|[;;]|[;;;]|[;;;;]|[;;;;;]|[\n]',line)
        NestVal[cnt-1,:3]=np.array([float(spltcont[1]),float(spltcont[2]), float(spltcont[3])])
        if cnt%NoNest==0:
            ThMatP1P4[4*(int(cnt/4)-1):4*int(cnt/4),:4]=nestposetotransfmat(NestVal[cnt-NoNest,:],NestVal[cnt-NoNest+1,:],NestVal[cnt-NoNest+2,:])      
print('Nest coordinates',NestVal,'Transformation matrices',ThMatP1P4)
transvec,toolposit=toolposit(ThMatP1P4[:4,:],ThMatP1P4[4:8,:],ThMatP1P4[8:12,:],ThMatP1P4[12:,:])
print('TransVector',transvec,'Tool Position',toolposit)
np.savetxt('TransToolPosit.csv',np.transpose(np.hstack((transvec,toolposit))))
terr=calctooldeferr(ThMatP1P4,transvec,toolposit)
print('Tool error',terr,'Mean error',np.mean(terr))
#Orientation def
fnamorient='ToolOrientDef.xyz'
OrientPosVal=np.zeros((3,3))
with open(fnamorient) as f:
    cnt=0
    for line in f:
        cnt+=1
        spltcont=re.split('[;]|[;;]|[;;;]|[;;;;]|[;;;;;]|[\n]',line)
        OrientPosVal[cnt-1,:3]=np.array([float(spltcont[1]),float(spltcont[2]), float(spltcont[3])])

xax=(OrientPosVal[1,:]-OrientPosVal[0,:])/np.linalg.norm(OrientPosVal[1,:]-OrientPosVal[0,:])#row 0 to row 1 is x axis
print('X axis vector',xax)
zax=(OrientPosVal[2,:]-OrientPosVal[0,:])/np.linalg.norm(OrientPosVal[2,:]-OrientPosVal[0,:])#row 0 to row 2 is z axis
print('z axis vector',zax)
yax=np.cross(zax,xax)/np.linalg.norm(np.cross(zax,xax))#X axis is assumed correct and zax is assumed to be in the XZ plane
zax=np.cross(xax,yax)/np.linalg.norm(np.cross(xax,yax))# actual z axis ensuring Orthonormality
RmatTractoFlan=np.transpose(np.vstack((np.vstack((xax,yax)),zax)))
#Considering pose 4
RmatFlantoNest=np.matmul(np.transpose(RmatTractoFlan),ThMatP1P4[:3,:3])
transFlangetoNest=-np.matmul(RmatFlantoNest,transvec)
ThFlantoNest=np.hstack((RmatFlantoNest,np.array([transFlangetoNest[0],transFlangetoNest[1],transFlangetoNest[2]])))
print('Flange to Nest coord',ThFlantoNest)
np.savetxt('TransFlangtoNest.csv',ThFlantoNest)
