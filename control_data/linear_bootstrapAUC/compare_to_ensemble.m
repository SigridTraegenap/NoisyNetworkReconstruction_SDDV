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
%where data are saved
sim_path = strcat(newdir,'/proof_of_concept/linear_sigma/');
save_string=sprintf('simulations/data_linear_T%d_dt%0.0e_Nnodes%d_Nincoming%d_sigma%d-%d',...
    T,delta_t, num_nodes, n_incoming, sigmamin, sigmamax);
save_string = strcat(sim_path, save_string);

Nnoise = 50;
Ntotal=10;

%% generate ensemble of connectivity matrices
Nens=500;

all_W = zeros(Nens, num_nodes, num_nodes);
for iens=1:Nens
    alpha=couplings(1);
    [betas,adjacency]=setup_NW(num_nodes, n_incoming, coupling, beta)
    all_W(iens,:,:)=adjancency;
end

%% loop over timeseries, reconstruct and compare to ensemble matrices (inner loop)
% flattening


all_AUCS = zeros(Nnoise, Ntotal, Nens);

for is=1:size(pre_s,2) 
    for irep=1:Ntotal
        load(sprintf(strcat(save_string, "_S%d_I%d.mat"), is, irep));  
         if sum(isnan(x_all(:)))>0
            all_AUCS(is, irep,:)=NaN;
         else
            %reconstruct
            [df_reconstr]=deriv_steadystate_reconstr(x_all)
            %loop over ensemble
            for iens=1:Nens
            [AUCs]=perfcurve_reconstr(df_reconstr, all_W(iens,:,:), num_nodes)
            all_AUCS(is, irep, iens)=AUCs
            end
         end
    end
end

all_AUCS = reshape(all_AUCS,Nnoise*Ntotal*Nens); 

save(sprintf('AUC_ensemble_linear_sigma%d-%d_dt%0.0e_res%0.0e',...
    sigmamin, sigmamax, delta_t, res), ...
    'all_AUCS')


function [AUCs]=perfcurve_reconstr(df_reconstr, adjacency, num_nodes)
bin_adjacency = adjacency;
bin_adjacency(bin_adjacency>0)=1;
bin_adjacency= reshape(bin_adjacency, num_nodes*num_nodes,1);
df_reconstr=reshape(df_reconstr, num_nodes*num_nodes,1);                                         
[~,~,~,AUCs] = perfcurve(bin_adjacency,...
                abs(df_reconstr), 1);
end


function [df_reconstr]=deriv_steadystate_reconstr(x_all)
[num_nodes, ~]=size(x_all);
df_reconstr=nan(num_nodes, num_nodes);
for inode2=1:num_nodes
   act_nextstep= x_all(inode2, 2:end);
   act_prestep = x_all(1:end,1:end-1);
   act_PINV = pinv(act_prestep); 
   df_reconstr(inode2,:)=mtimes(act_nextstep,act_PINV);    
end
df_reconstr(logical(eye(size(df_reconstr)))) = 0.;
end

%nextstep