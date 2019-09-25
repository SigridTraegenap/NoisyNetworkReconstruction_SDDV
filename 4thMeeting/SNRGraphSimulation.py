"""
Created on Fri Aug 24 08:52:00 2018

@author: Jose Casadiego
"""


#sys.path.insert(0, 'C:\\Users\\casadieg\\Documents\\Projects\\Simulation_Suite\\Python\\')
#sys.path.insert(0, '/home/beate/Dokumente/00. UNI/0.Stipendium/Kolleg/network_reconstruction/4.Treffen/Data_Jose_SNR/Noise_driven/GoodDatawithJandA/')
import numpy as np
import _pickle as pickle
#import Graphics as gr
import os
from sklearn import metrics
import sys  
import itertools as it
#clear()

J=np.loadtxt("/home/beate/Dokumente/00. UNI/0.Stipendium/Kolleg/network_reconstruction/4.Treffen/Data_Jose_SNR/Noise_driven/GoodDatawithJandA/connectivity.dat")
b=np.eye(np.size(J))
#Ad=np.copy(J)
#Ad[Ad!=0]=1
#Ad=Ad+np.eye(len(Ad))
#data=pickle.load( open("Data/data.p", "rb" ) )

files=os.listdir('/home/beate/Dokumente/00. UNI/0.Stipendium/Kolleg/network_reconstruction/4.Treffen/Data_Jose_SNR/Noise_driven/GoodDatawithJandA/Data/.')
#print(files)
data=0
for i in range(1,len(files)):
    k=data
    #print('k',np.isnan(k))
    with open("/home/beate/Dokumente/00. UNI/0.Stipendium/Kolleg/network_reconstruction/4.Treffen/Data_Jose_SNR/Noise_driven/GoodDatawithJandA/Data/%s"%files[i], "rb" ) as fd:
        data=pickle.load(fd, fix_imports=True, encoding="latin1")
        #data = np.asarray(data)
        #np.nan_to_num(data,copy=True)
        #print(data)
        #data = np.ndarray.tolist(data)
        A = np.asarray(data)
 
        #A=data
        exp=np.zeros(20)
        for l in range(0,20): #if l are the time nodes
            for k in range(0,20):
                a=np.reshape(A[l,k,:],20000)
                x=np.asarray(a)
                x3=a**3
                xexphelp=x3*np.transpose(x)
                matr=np.reshape(xexphelp,np.size(xexphelp))
                exp[l]=np.mean(matr)
                print(exp)
        
            '''      
        
        
        
        for l in range(0,20): #if l are the time nodes
            a=np.reshape(A[l,:,:],20000)
            u=np.ndarray.tolist(a)
            d = {x:u.count(x) for x in u}    
            all_keys=[]
            all_items=[]
            for i,j in d.items():
                all_keys.append(i)
                all_items.append(j)
            print(type(all_keys),type(all_items))
            #occ, val = d.keys(), d.values()
            #b1=list(occ) #values that are counted
            #b2=list(val) #how often they come
            b1=np.asarray(all_keys)
            b2=np.asarray(all_items)
            b2=b2/20000 #probabilities
            #now get expectation value of x**3 * x^T
            exp[l]=np.sum(b1*b2)
            
        #print(data.shape)
        
        #print(data)
        
        
    #data=pickle.load( open(  )
        
        #print ('Working with file: %s'%files[i])
       # print(np.mean(data))
       

    
    T=len(data)
    N,M=data[0].shape
    AUCs=np.zeros((N,T))
    
    for t in range(T):
        print "Ensemble: %d/%d"%(t,T)
        x=np.zeros((N,(M-1)*T))
        y=np.zeros((N,(M-1)*T))
        for r in range(t):
            x[:,r*(M-1):(r+1)*(M-1)]=data[r][:,:-1]
            y[:,r*(M-1):(r+1)*(M-1)]=data[r][:,1:]
        
        try:
            L=np.dot(y,np.linalg.pinv(x))
        
        except np.linalg.LinAlgError:
            L=np.zeros((N,N))
            
        
        for n in range(N):
            fpr, tpr, thresholds = metrics.roc_curve(Ad[n,:],np.fabs(L[n,:]),pos_label=1)
            AUCs[n,t]=metrics.auc(fpr, tpr)
        
    np.savetxt("Results/reconstruction_%s.dat"%files[i][-5:-2],AUCs,delimiter="\t",fmt="%1.6f")
    '''