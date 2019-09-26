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
sigma = 10^-5;
sigma_add=1/6*10^-3;

%where to save data
save_string=sprintf('simulations/data_linear_T%d_dt%0.0e_Nnodes%d_Nincoming%d_sigma%0.0e',...
    T,delta_t, num_nodes, n_incoming, sigma);
Ntotal=50;

%% simulate N timeseries without noise
% simulate N timeseries with noise using same connectivity
alpha=couplings(1);
Ntotal=50;
for isim=1:Ntotal
    [x_all, x_tau_all, dt_x_all, dt_x_all_v2, adjacency, betas]  = simulate_lineardecay_tseries(...
                    sigma,num_simulations,alpha,num_nodes,...
                    n_incoming, T, delta_t, res,initial);
    save(sprintf(strcat(save_string, "precisiontest_I%d.mat"), isim));

end

%% without rounding or noise

all_AUCS=zeros(Ntotal,1);
for isim=1:Ntotal
    load(sprintf(strcat(save_string, "precisiontest_I%d.mat"), isim));
    [AUCs]=reconstruction_2p_approx(x_tau_all,...
                                   dt_x_all, adjacency,num_nodes);
    all_AUCS(isim)=AUCs;
end

save(sprintf("simulations/summary_2palg_precisionfull_S%0.0e.mat", sigma),  'all_AUCS');

%% with rounding

all_AUCS=zeros(Ntotal,1);
for isim=1:Ntotal
    load(sprintf(strcat(save_string, "precisiontest_I%d.mat"), isim));
    [AUCs]=reconstruction_2p_approx(round(x_tau_all,3),...
                                   round(dt_x_all,3), adjacency,num_nodes);
    all_AUCS(isim)=AUCs;
end

save(sprintf("simulations/summary_2palg_precisionround_S%0.0e.mat", sigma),  'all_AUCS');

%% with additional noise

all_AUCS=zeros(Ntotal,1);
for isim=1:Ntotal
    load(sprintf(strcat(save_string, "precisiontest_I%d.mat"), isim));
    [AUCs]=reconstruction_2p_approx(x_tau_all + sigma_add*randn(size(x_tau_all)),... 
                                   dt_x_all+sigma_add*randn(size(dt_x_all)),...
                                   adjacency,num_nodes);
    all_AUCS(isim)=AUCs;
end

save(sprintf("simulations/summary_2palg_precisionaddnoise_S%0.0e.mat", sigma),  'all_AUCS');

%% bar/scatter plot

% comment: if I don't average AUC across nodes for one reconstruction, I
% get very low and very high AUC (altlugh it should be around 0.5??) 
%I now try and average across nodes
%previous: reshape(all_AUCS, Ntotal*num_nodes,1);

%load wo loss of precision
load(sprintf("simulations/summary_2palg_precisionfull_S%0.0e.mat", sigma))
auc_wonoise = all_AUCS;
%load with rounding
load(sprintf("simulations/summary_2palg_precisionround_S%0.0e.mat", sigma))
auc_round = all_AUCS;
%load with additional noise
load(sprintf("simulations/summary_2palg_precisionaddnoise_S%0.0e.mat", sigma))
auc_addnoise = all_AUCS;

xdata1 = ones(size(auc_wonoise)) + randn(size(auc_wonoise))*0.1;
xdata2 = ones(size(auc_wonoise))*2 + randn(size(auc_wonoise))*0.1;
xdata3 = ones(size(auc_wonoise))*3 + randn(size(auc_wonoise))*0.1;

figure('Name', 'BarScatter_AU')
scatter(xdata1, auc_wonoise)
hold on
scatter(1, mean(auc_wonoise), "ks", 'markerfacecolor', 'k')
hold on
scatter(xdata2, auc_round)
hold on
scatter(2, mean(auc_round), "ks", 'markerfacecolor', 'k')
hold on
scatter(xdata3, auc_addnoise)
hold on
scatter(3, mean(auc_addnoise), "ks", 'markerfacecolor', 'k')
ylabel("AUC")
%xticks([1,2]);
%xticklabels({"w\\o noise", "w\ noise"})