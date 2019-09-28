%setup
clear;
close all;

%add code database
mydir  = pwd;
idcs   = strfind(mydir,filesep);
newdir = mydir(1:idcs(end-1)-1);
addpath(newdir + "/reconstruction_simulation_networks/") 

%define parameters
num_nodes=30; %number_nodes
n_incoming=5; %number incoming nodes (sparsity)
delta_t=0.001;  %fine simulation tscale
res=0.01;  %coarse reconstruction t-scale
initial="random"; %specify initial conditions
T=100;  %how many timesteps are evaluated
couplings=[1];  %coupling strength 
num_simulations=1;  %how many simulations to perform 

%noise
sigmamin = -5;
sigmamax = 3;
%precision

sigma_add=1/6*10^-3;

%where to save data


Nnoise = 50;
Ntotal=10;


simulate=false;
delta_t=0.0001;
res=0.01;

save_string=sprintf('simulations/data_cubic_T%d_dt%0.0e_Nnodes%d_Nincoming%d_sigma%d-%d',...
    T,delta_t, num_nodes, n_incoming, sigmamin, sigmamax);

pre_s=linspace(sigmamin,sigmamax,Nnoise);
pre_s= power(10,pre_s);  %log-spacing of noise
alpha=couplings(1);
for is=1:3
    sigma=pre_s(is);   
    for irep=1:Ntotal
        [x_all, x_tau_all, dt_x_all, dt_x_all_v2, adjacency, betas]=simulate_cubicdecay_tseries(...
            sigma,num_simulations,alpha,num_nodes,...
            n_incoming, T, delta_t, res,initial);
        save(sprintf(strcat(save_string, "_S%d_I%d.mat"), is, irep));    
    end
end




