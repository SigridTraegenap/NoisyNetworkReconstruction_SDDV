#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 24 14:10:01 2019

@author: beate
"""

from sklearn import metrics
import numpy as np
import pylab as plt
import _pickle as pickle
import scipy.io
import os
from reconstruction_radio_dt_res import simulate
'''
here I calculate the meanAUC (and therefore how good my reconstruction is) for alphas between 10^-5 and 10^2 for different corase to fine scale ratios.
Therefore in the end i have a 3dim matrix 
R ... the save vector for the AUC values
c ... the factor in the reconstruction file of the model function

'''
plt.close() #he closes the plots from before
#creating noise vector
Numb=80
sigma=np.linspace(-5,2,Numb)
sigma=10**sigma

#creating radio vector
resvec=np.linspace(0.001,10,10)

#for saving stuff
r=np.random.rand(1) #this is just so the files don't overwrite each other because I want to plot more than one random data file
#cond='steady state' #one of the two initial conditions
cond='random state'
b='random'
#setting simulation parameters
c=1 #power fe cubic or linear
N=30 #number of units. = Knoten   
o=2 #time interval for simulation [0,M].
T=1 #time series

Z=np.zeros((2,Numb+1,np.size(resvec))) #save matrix

for l in range(len(resvec)):
    R=np.zeros(np.size(sigma)) #array to save AUCs in 
    for i in range(len(sigma)):
    
        print(resvec[l])
        simulate(N,o,sigma[i],T,c,resvec[l])
        #Code from reconstruct_expectations.py
    
        J=np.zeros((N,N))
        data=np.around(pickle.load( open("Data/data.p", "rb" ) ),4)#around rounds the datapoints to 4 signifitcant 
        J=np.loadtxt("Data/connectivity.dat")
        Ad=np.copy(J)
        Ad[Ad!=0]=1
        
        Ad=Ad+np.eye(len(Ad))
        
        T=len(data)
        N,M=data[0].shape
        
        #T=10 #Verschiedene Timeserien
        
        x=np.zeros((N,(M-1)*T)) #Zeitserie jetzt
        y=np.zeros((N,(M-1)*T)) #Zeitserie zukunft
        for t in range(T):
            x[:,t*(M-1):(t+1)*(M-1)]=data[t][:,:-1]
            y[:,t*(M-1):(t+1)*(M-1)]=data[t][:,1:]
       # print(x)    
        L=np.dot(y,np.linalg.pinv(x)) #compute pseudo inverse of x
        print('AUC step')
        AUCs=np.zeros((N,1))
        for n in range(N):
            fpr, tpr, thresholds = metrics.roc_curve(Ad[n,:],np.fabs(L[n,:]),pos_label=1)
            AUCs[n,0]=metrics.auc(fpr, tpr)
    
        R[i]=np.mean(AUCs)
        #print(R)
     #3 columns and 80!! is the number adaptable to the number of sigmas
    Z[0,:]=sigma
    Z[1,:]=R 
    Z[2,0]=resvec[l]
    np.savetxt("Data/AUCvalues/AUC_c{}_initalconditions{}_beta{}_M{}_N{}_T{}_res{}.dat".format(c,cond,b,M,N,T,resvec[l]),Z,delimiter="\t",fmt="%1.4f")          
#print(R)    

    