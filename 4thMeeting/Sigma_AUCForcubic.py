from sklearn import metrics
import numpy as np
import pylab as plt
import _pickle as pickle
import scipy.io
import os
from resconstruction_cubic import simulate
'''
here I calculate the meanAUC (and therefore how good my reconstruction is) for alphas between 10^-5 and 10^2.
R ... the save vector for the AUC values
c ... the factor in the reconstruction file of the model function
in the end it plots it
'''
plt.close() #he closes the plots from before
#worked for c=3 and c=1, but for c=5 doesn't converge so we take another one
Numb=80
sigma=np.linspace(-5,2,Numb)
sigma=10**sigma
'''
sigma=np.zeros(90)
Numb=90 #because Z is defined by this size
Numbb=80
sigma1=np.linspace(-5,1,Numbb)
sigma1=10**sigma1
sigma2=np.linspace(11,30,10) #it seems that the network is driven too far away to converge again if i put in values above 30
sigma[0:80]=sigma1
sigma[80:90]=sigma2
'''
#for saving stuff
r=np.random.rand(1) #this is just so the files don't overwrite each other because I want to plot more than one random data file
#cond='steady state' #one of the two initial conditions
cond='random state'
b='random'
#setting simulation parameters
c=3 #power fe cubic or linear
N=30 #number of units. = Knoten   
o=10 #time interval for simulation [0,M].
T=1 #time series

R=np.zeros(np.size(sigma)) #array to save AUCs in 
for i in range(len(sigma)):

    print('sigma',sigma[i])
    simulate(N,o,sigma[i],T,c)
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
Z=np.zeros((2,Numb)) #2 columns and 80!! is the number adaptable to the number of sigmas
Z[0,:]=sigma
Z[1,:]=R 
np.savetxt("Data/AUCvalues/AUC_c{}_initalconditions{}_beta{}_M{}_N{}_T{}_random{}.dat".format(c,cond,b,M,N,T,r),Z,delimiter="\t",fmt="%1.4f")          
#print(R)    
#
    
plt.figure()
plt.plot(sigma,R)

plt.ylabel('AUC')
plt.xlabel('sigma')
plt.title('Plot wuhu')
plt.show()