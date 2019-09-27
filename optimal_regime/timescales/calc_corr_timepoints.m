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

%%
%where data are saved
sim_path = strcat(newdir,'/proof_of_concept/linear_sigma/');
save_string=sprintf('simulations/data_linear_T%d_dt%0.0e_Nnodes%d_Nincoming%d_sigma%d-%d',...
    T,delta_t, num_nodes, n_incoming, sigmamin, sigmamax);
save_stringGEN = strcat(sim_path, save_string);

Nnoise = 50;
Ntotal=10;

%% x_all, addnoise

for is=1:Nnoise 
    for irep=1:Ntotal
        load(sprintf(strcat(save_stringGEN, "_S%d_I%d.mat"), is, irep), 'x_all');  
        x_all = x_all + sigma_add*randn(size(x_all));
         if sum(isnan(x_all(:)))>0
            all_AC(is, irep,:)=NaN;
         else
            parfor inode=1:num_nodes                
                all_AC(is, irep,inode)=custom_autocorr(squeeze(x_all(inode,:)), 1);
            end           
         end
    end
end

save(sprintf('summary_linear_autocorr_xall_sigma%d-%d_dt%0.0e_res%0.0e',...
    sigmamin, sigmamax, delta_t, res), 'all_AC')


%%
%where data are saved

sim_path = strcat(newdir,'/optimal_regime/vary_sigma_cubic/');
save_string=sprintf('simulations/data_cubic_T%d_dt%0.0e_Nnodes%d_Nincoming%d_sigma%d-%d',...
    T,delta_t, num_nodes, n_incoming, sigmamin, sigmamax);
save_stringGEN = strcat(sim_path, save_string);

Nnoise = 50;
Ntotal=10;

%% x_all, addnoise
for is=1:Nnoise 
    for irep=1:Ntotal
        load(sprintf(strcat(save_stringGEN, "_S%d_I%d.mat"), is, irep), 'x_all');  
        x_all = x_all + sigma_add*randn(size(x_all));
         if sum(isnan(x_all(:)))>0
            all_AC(is, irep,:)=NaN;
         else
            parfor inode=1:num_nodes
                all_AC(is, irep,inode)=custom_autocorr(squeeze(x_all(inode,:)), 1);
            end           
         end
    end
end

save(sprintf('summary_cubic_autocorr_xall_sigma%d-%d_dt%0.0e_res%0.0e',...
    sigmamin, sigmamax, delta_t, res), 'all_AC')

%%
function ac=custom_autocorr(y, lag)
y1=y(lag+1:end);
y2=y(1:end-lag);
ac=corr(y1.',y2.');
end
