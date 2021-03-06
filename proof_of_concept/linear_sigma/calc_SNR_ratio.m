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
save_string=sprintf('simulations/data_linear_T%d_dt%0.0e_Nnodes%d_Nincoming%d_sigma%d-%d',...
    T,delta_t, num_nodes, n_incoming, sigmamin, sigmamax);

Nnoise = 50;
Ntotal=10;
pre_s=linspace(sigmamin,sigmamax,Nnoise);
pre_s= power(10,pre_s);  %log-spacing of noise

%% collect SNR, rounding
all_sd = zeros(Nnoise, Ntotal, num_nodes);

for is=1:Nnoise 
    var=pre_s(is)^2;
    for irep=1:Ntotal
        load(sprintf(strcat(save_string, "_S%d_I%d.mat"), is, irep));  
        x_all_noise = round(x_all(:,int16(end/2):end), 3);
        sds = (std(x_all_noise,0,2).^2 + mean(x_all_noise, 2).^2)/var;
        all_sd(is, irep,:)=sds;
    end
end

save(sprintf('simulations/summary_std_sigma_round_sigma%d-%d', sigmamin, sigmamax), ...
    'all_sd')

%% collect SNR, add noise
all_sd = zeros(Nnoise, Ntotal, num_nodes);

for is=1:Nnoise 
    var=pre_s(is)^2;
    for irep=1:Ntotal
        load(sprintf(strcat(save_string, "_S%d_I%d.mat"), is, irep));  
        x_all_noise = x_all(:,int16(end/2):end);
        x_all_noise = x_all_noise + sigma_add*randn(size(x_all_noise));
        sds = (std(x_all_noise,0,2).^2 + mean(x_all_noise, 2).^2)/var;
        all_sd(is, irep,:)=sds;
    end
end

save(sprintf('simulations/summary_std_sigma_addnoise_sigma%d-%d', sigmamin, sigmamax), ...
    'all_sd')

%% collect SNR, full precision
all_sd = zeros(Nnoise, Ntotal, num_nodes);

for is=1:Nnoise 
    var=pre_s(is)^2;
    for irep=1:Ntotal
        load(sprintf(strcat(save_string, "_S%d_I%d.mat"), is, irep));  
        x_all_noise = x_all(:,int16(end/2):end);        
        sds = (std(x_all_noise,0,2).^2 + mean(x_all_noise, 2).^2)/var;
        all_sd(is, irep,:)=sds;
    end
end

save(sprintf('simulations/summary_std_sigmafullprec_sigma%d-%d', sigmamin, sigmamax), ...
    'all_sd')