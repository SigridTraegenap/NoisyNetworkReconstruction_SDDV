# -*- coding: utf-8 -*-
"""
Created on Fri Jun 29 17:26:24 2018

@author: Jose Casadiego
"""

import numpy as np
import pylab as pl
import _pickle as pickle
import scipy.io
import os


def topology(N,ni):
    """
    topology(N,ni) generates connectivity matrices for network simulation.
    
    Parameters:
    ------------------
    N:  number of units.
    ni: number of incoming connections per unit.
    
    Ouput:
    ------------------
    J:  Weigthed adjacency matrix with (at most) ni incoming links per unit. 
        The ouput is also saved in "Data/connectivity.dat".
    """
    
    J=np.zeros((N,N))   
    while sorted(np.linalg.eigvals(J))[2] == 0: #we define that it has to be strongly connected
        for i in range(N):
            f=np.random.permutation(N)[:ni] #random.üermutations creates N random nunbers in a vector, [:ni] takes the first ni numbers of this vector
           #J[i,f]=0.2+0.8*np.random.rand(1,len(f)) #len = returning number of elements in "f
            J[i,f]=np.random.randint(2)*0.5 # creates randomly either 0 or 1
        np.fill_diagonal(J,0)
        np.savetxt("Data/connectivity.dat",J,delimiter="\t",fmt="%1.4f")
       # a=np.linalg.eigvals(J)
    return J

def beta(N):
    """
    beta(N) generates the degradation rate for each unit in the network.
    
    Parameters:
    ------------------
    N:  number of units.
    
    Ouput:
    ------------------
    b:  Vector of random values uniformly distributed in [0.2,1.0]. 
    b=np.zeros(N)
    for i in range(0,N):
        b[i]=1
    """
    
        
    b=1
    #b=0.2+0.8*np.random.rand(N)
    
    return b

def simulate(N,M,sigma,T,J,b):
    """
    simulate(N,M,sigma,T) generate noisy simulations of networks of dynamical 
    systems with a fixed number of incoming connections per unit.
    
    Parameters:
    ------------------
    N:      number of units. = Knoten
    M:      time interval for simulation [0,M].
    sigma:  strength of noisy signal.
    T:      number of time series to be simulated from different initial conditions
            
    Ouput:
    ------------------
    Data/dynamics.p:    pickled file containing the dynamics for all
                        the simulated conditions.
    Figure:             plot demonstrating how the dynamics of two different
                        units look in different sampling scales
                    
    Example:
    ------------------
    simulate(30,10,0.6,4) simulates a network of 30 units having 5 incoming 
    links per unit in the interval [0,10] with a noise strength of 0.6 using
    4 different initial conditions (thus leading to 4 different time series).
    
    """
    
    pl.close("all")
    pl.style.use("ggplot")
    
    #Defining network parameters
    #J=J
    #b=b
    
    #Defining simulation paramerters
    res=0.01             #Coarse scale 
    dt=0.001            #Original scale
    sqrtdt=np.sqrt(dt)
    
    
    tspan=np.arange(0,M,res)    #Sampling interval in coarse scale
    tspan2=np.arange(0,M,dt)    #Sampling interval in original scale
    
    L=len(tspan)
    L2=len(tspan2)
    
    ind=range(0,L2,int(res/dt)) #Sampling indices
    
    data=[]
    
    print ("Simulating:")
    
    for t in range(T):
        
        print ("Time series: %d/%d"%(t+1,T))
        
        Y=np.zeros((N,L))       #Dynamics in coarse scale
        y=np.zeros((N,L2))      #Dynamics in original scale
        
        #y[:,0]=-2.0 + 4.0*np.random.rand(N)     #Random initial conditions
        y[:,0]=0
        #Simulation of noisy dynamics
        for m in range(L2-1):
            
            S=np.random.rand(N,1)
            S[S>=0.5]=1
            S[S<0.5]=-1
        
            W=sqrtdt*np.random.randn(N,1)
            K1=dt*model(m,y[:,m],J,b)+sigma*(W-sqrtdt*S)
            K2=dt*model(m,y[:,m]+K1.reshape(N),J,b)+sigma*(W+sqrtdt*S)
            
            y[:,m+1]=y[:,m]+0.5*(K1+K2).reshape(N)
            
        Y=y[:,ind]
        data.append(Y)
    
    #Plotting the dynamics of two units in the two scales
    '''
    f, (ax1, ax2) = pl.subplots(2)
    ax1.plot(y[5:7,:].T)
    ax1.set_title("Sampling the dynamics in two different scales\nOriginal scale")
    ax1.set_ylabel(r"$x_{i}(t)$")
    ax2.plot(Y[5:7,:].T)
    ax2.set_title("Coarser scale")
    ax2.set_xlabel("Iteration")
    ax2.set_ylabel(r"$x_{i}(t)$")
    pl.tight_layout() '''
    
    #Saving data files
    directory=os.getcwd()
    pickle.dump(data,open("%s/Data/data.p"%directory,"wb"))
    scipy.io.savemat("%s/Data/data.mat"%directory, mdict={'data': data})
    print ("Simulations completed.\nOutputs saved to:\n 'Data/data.p' (Pickled data)\n 'Data/connectivity.dat' (Numpy array)")

def model(t,y,J,b):
    """
    model(t,y,J,b) computes the rate of change of network dynamical systems
    described by:
        
        \dot{x}_i(t)=-b_i*x_i(t)+sum_{i=1}^{N}J_{ij}*(x_j(t)-x_i(t))
    
    Parameters:
    ------------------
    t:  time step.
    y:  vector containing network state at time step t.
    J:  connectivity matrix.
    b:  vector containing degradation rates for every unit.
            
    Ouput:
    ------------------
    dydt:   vector containing the rate of change of every unit in the network
    
    """
   
    N=len(J)
    
    dydt=np.zeros((N,1))  

    for i in range(N):
        suma=0.0
        for j in range(N):
            suma=J[i,j]*(y[j]-y[i])+suma
        #dydt[i]=-b[i]*y[i]+suma
        dydt[i]=-b*y[i]**3+suma
          
    return dydt


#Running an example
if __name__=="__main__": #only do simulate if the program is called explicitly 
    J=topology(40,5)
    simulate(40,10,0,2,J,1) 