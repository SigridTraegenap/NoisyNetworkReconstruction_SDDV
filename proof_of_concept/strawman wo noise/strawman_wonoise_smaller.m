%setup
clear;
close all;

%add code database
mydir  = pwd;
idcs   = strfind(mydir,filesep);
newdir = mydir(1:idcs(end-1)-1);
addpath(newdir + "/reconstruction_simulation_networks/") 

%define parameters
num_nodes=10; %number_nodes
n_incoming=2; %number incoming nodes (sparsity)
delta_t=0.001;  %fine simulation tscale
res=0.01;  %coarse reconstruction t-scale
initial="random"; %specify initial conditions
T=100;  %how many timesteps are evaluated
couplings=[1];  %coupling strength 
num_simulations=1;  %how many simulations to perform 

%noise
sigma = 10^-1;
Ntotal=50;

%where to save data
save_string=sprintf('simulations/data_linear_T%d_dt%0.0e_Nnodes%d_Nincoming%d_sigma%0.0e',...
    T,delta_t, num_nodes, n_incoming, sigma);
%% simulate N timeseries without noise
% simulate N timeseries with noise using same connectivity
alpha=couplings(1);

for isim=1:Ntotal
    rng(isim)
    [x_all, x_tau_all, dt_x_all, dt_x_all_v2, adjacency, betas]  = simulate_lineardecay_tseries(...
                    0,num_simulations,alpha,num_nodes,...
                    n_incoming, T, delta_t, res,initial);
    save(sprintf(strcat(save_string, "wo_noise_I%d.mat"), isim));

    rng(isim)
    [x_all, x_tau_all, dt_x_all, dt_x_all_v2, adjacency, betas]  = simulate_lineardecay_tseries(...
                    sigma,num_simulations,alpha,num_nodes,...
                    n_incoming, T, delta_t, res,initial, betas, adjacency);
    save(sprintf(strcat(save_string, "w_noise_I%d.mat"), isim));
end


%% collect results wo noise
all_AUCS=zeros(Ntotal,1);
all_dfs = zeros(Ntotal, num_nodes, num_nodes);
for isim=1:Ntotal
    load(sprintf(strcat(save_string, "wo_noise_I%d.mat"), isim));
    [df_reconstr]=deriv_steadystate_reconstr(round(x_tau_all,3),...
                                             round(dt_x_all,3) );
    all_dfs(isim,:,:)=df_reconstr;
    [AUCs]=reconstruction_2p_approx(round(x_tau_all,3), ...
                    round(dt_x_all,3), adjacency,num_nodes);
    all_AUCS(isim)=AUCs;
end

save("simulations/summary_smaller_wo_noise_2palg.mat", 'all_dfs', 'all_AUCS');


%% collect results w noise
all_AUCS=zeros(Ntotal,1);
all_dfs = zeros(Ntotal, num_nodes, num_nodes);
for isim=1:Ntotal
    load(sprintf(strcat(save_string, "w_noise_I%d.mat"), isim));
    [df_reconstr]=deriv_steadystate_reconstr(round(x_tau_all,3),...
                                             round(dt_x_all,3) );
    all_dfs(isim,:,:)=df_reconstr;
    [AUCs]=reconstruction_2p_approx(round(x_tau_all,3), ...
                    round(dt_x_all,3), adjacency,num_nodes);
    all_AUCS(isim)=AUCs;
end
disp(mean(all_AUCS(:)));
save("simulations/summary_smaller_w_noise_2palg.mat", 'all_dfs', 'all_AUCS');

%% collect results wo noise - nextstep
all_AUCS=zeros(Ntotal,1);
all_dfs = zeros(Ntotal, num_nodes, num_nodes);
for isim=1:Ntotal
    load(sprintf(strcat(save_string, "wo_noise_I%d.mat"), isim));
    [df_reconstr]=Nextstep_deriv_steadystate_reconstr(x_all );
    all_dfs(isim,:,:)=df_reconstr;
    [AUCs]=reconstruction_nextstep(x_all, ...
                                    adjacency,num_nodes);
    all_AUCS(isim)=AUCs;
end

save("simulations/summary_smaller_wo_noise_nextstepalg.mat", 'all_dfs', 'all_AUCS');


%% collect results w noise - nextstep
all_AUCS=zeros(Ntotal,1);
all_dfs = zeros(Ntotal, num_nodes, num_nodes);
for isim=1:Ntotal
    load(sprintf(strcat(save_string, "w_noise_I%d.mat"), isim));
    [df_reconstr]=Nextstep_deriv_steadystate_reconstr(x_all );
    all_dfs(isim,:,:)=df_reconstr;
    [AUCs]=reconstruction_nextstep(x_all, ...
                                    adjacency,num_nodes);
    all_AUCS(isim)=AUCs;
end
disp(mean(all_AUCS(:)));
save("simulations/summary_smaller_w_noise_nextstepalg.mat", 'all_dfs', 'all_AUCS');

%% plot connectivity  matrix (single) and reconstructed matrix

isim=1;
%load original
load(sprintf(strcat(save_string, "wo_noise_I%d.mat"), isim));
adjacency_true = adjacency;

%load wo noise reconstructed
load("simulations/summary_smaller_wo_noise_nextstepalg.mat")
adjacency_wonoise = squeeze(all_dfs(isim,:,:));
%load w noise reconstructed
load("simulations/summary_smaller_w_noise_nextstepalg.mat")
adjacency_wnoise = squeeze(all_dfs(isim,:,:));

min1 = min(adjacency_true(:));
min2 = min(adjacency_wonoise(:));
min3 = min(adjacency_wnoise(:));
vmin =min(min(min1, min2), min3);

max1 = max(adjacency_true(:));
max2 = max(adjacency_wonoise(:));
max3 = max(adjacency_wnoise(:));
vmax =max(max(max1, max2), max3);

clims = [0,max(abs(vmax), abs(vmin))]; 

figure('Name', 'Reconstructed connectivity matrix') 
subplot(131)
imagesc(adjacency_true)
subplot(132)
imagesc(abs(adjacency_wonoise))
subplot(133)
imagesc(abs(adjacency_wnoise))
colorbar()


figure(2)
scatter(reshape(adjacency_true, num_nodes.^2,1), reshape(adjacency_wnoise,num_nodes.^2,1))
%% bar/scatter plot

% %load wo noise reconstructed
% load("simulations/summary_wo_noise_2palg.mat")
% auc_wonoise = reshape(all_AUCS, Ntotal*num_nodes,1);
% %load w noise reconstructed
% load("simulations/summary_w_noise_2palg.mat")
% auc_wnoise = reshape(all_AUCS, Ntotal*num_nodes,1);
% 
% xdata1 = ones(Ntotal*num_nodes,1) + randn(Ntotal*num_nodes,1)*0.1;
% xdata2 = ones(Ntotal*num_nodes,1)*2 + randn(Ntotal*num_nodes,1)*0.1;


%load wo noise reconstructed
load("simulations/summary_wo_noise_2palg.mat")
auc_wonoise = all_AUCS;
%load w noise reconstructed
load("simulations/summary_w_noise_2palg.mat")
auc_wnoise = all_AUCS;

xdata1 = ones(size(auc_wonoise)) + randn(size(auc_wonoise))*0.1;
xdata2 = ones(size(auc_wonoise))*2 + randn(size(auc_wonoise))*0.1;

figure('Name', 'BarScatter_AU')
scatter(xdata1, auc_wonoise)
hold on
scatter(1, mean(auc_wonoise), "ks", 'markerfacecolor', 'k')
hold on
scatter(xdata2, auc_wnoise)
hold on
scatter(2, mean(auc_wnoise), "ks", 'markerfacecolor', 'k')
ylabel("AUC")
xticks([1,2]);
xticklabels({"w\\o noise", "w\ noise"})
