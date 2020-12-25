# Laser Tracker Pose Data for Kinematic Identification
The data used as part of the article with the following details: R. A. Boby, A. Klimchik, Combination of Geometric and Parametric Approaches for Kinematic Identification of an Industrial Robot, Robotics and Computer Integrated Manufacturing (under review). Data includes the end-effector pose data measured using laser tracker (N1N2N3Jval.csv), identified values of Tool to Flange transformation matrix (TransFlangtoNest0.csv) and Base to Robot coordinate frames (TransTractoWorld0.csv). Details of each files are as follows:

1. End-effector pose data (N1N2N3Jval.csv)

The first 9 columns of N1N2N3Jval.csv represent the position data of the 3 SMR mounted on the end-effector of the robot measured using a laser tracker. One can use 3 postion data to compute the position and orientation of the coordinate frame defined to align with the SMR system. 
The last 9 columns represent the joint values from the robot controller.

2. Identified values of Tool to Flange transformation matrix (TransFlangtoNest0.csv)

The 3 x 4 matrix representing the transformation from Flange (F in article) to Tool coordinate systems (T in article).

3. Base to Robot coordinate frames (TransTractoWorld0.csv)

The 3 x 4 matrix representing the transformation from Tracker (B) to robot coordinate systems (R).

Python codes do the following:

1) NesttoFlangeV2_1.py: Helps in estimation of the tranfromation between tool and flanfe coordinate frames.

2) 
