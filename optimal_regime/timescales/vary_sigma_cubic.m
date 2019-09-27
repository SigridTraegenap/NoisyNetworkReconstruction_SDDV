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
save_string=sprintf('simulations/data_cubic_T%d_dt%0.0e_Nnodes%d_Nincoming%d_sigma%d-%d',...
    T,delta_t, num_nodes, n_incoming, sigmamin, sigmamax);

Nnoise = 50;
Ntotal=10;

pre_s=linspace(sigmamin,sigmamax,Nnoise);
pre_s= power(10,pre_s);  %log-spacing of noise
%% simulate timeseries


alpha=couplings(1);
for is=1:size(pre_s,2)
    sigma=pre_s(is);   
    for irep=1:Ntotal
        [x_all, x_tau_all, dt_x_all, dt_x_all_v2, adjacency, betas]=simulate_cubicdecay_tseries(...
            sigma,num_simulations,alpha,num_nodes,...
            n_incoming, T, delta_t, res,initial);
        save(sprintf(strcat(save_string, "_S%d_I%d.mat"), is, irep));    
    end
end

%% collect results for 2p algorithms, rounding
all_AUCS = zeros(Nnoise, Ntotal);

for is=1:size(pre_s,2)
 
    for irep=1:Ntotal
        load(sprintf(strcat(save_string, "_S%d_I%d.mat"), is, irep));  
         if sum(isnan(x_all(:)))>0
            all_AUCS(is, irep)=NaN;
         else
            [AUCs]=reconstruction_2p_approx(round(x_tau_all,3), ...
                        round(dt_x_all,3), adjacency,num_nodes);
            all_AUCS(is, irep)=AUCs;
         end
    end
end

save(sprintf('simulations/summary_2p_round_sigma%d-%d_dt%0.0e_res%0.0e',...
    sigmamin, sigmamax, delta_t, res), ...
    'all_AUCS')
%% collect results for 2p alg, add noise
all_AUCS = zeros(Nnoise, Ntotal);

for is=1:Nnoise
    for irep=1:Ntotal
        load(sprintf(strcat(save_string, "_S%d_I%d.mat"), is, irep));  
        if sum(isnan(x_all(:)))>0
            all_AUCS(is, irep)=NaN;
         else
        [AUCs]=reconstruction_2p_approx(x_tau_all + sigma_add*randn(size(x_tau_all)), ...
                    dt_x_all + sigma_add*randn(size(dt_x_all)), ...
                    adjacency,num_nodes);
        all_AUCS(is, irep)=AUCs;
        end
    end
end

save(sprintf('simulations/summary_2p_addnoise_sigma%d-%d_dt%0.0e_res%0.0e',...
    sigmamin, sigmamax, delta_t, res), ...
    'all_AUCS')

%% collect results for 3p algorithms, rounding
all_AUCS = zeros(Nnoise, Ntotal);

for is=1:size(pre_s,2)
 
    for irep=1:Ntotal
        load(sprintf(strcat(save_string, "_S%d_I%d.mat"), is, irep));    
        if sum(isnan(x_all(:)))>0
            all_AUCS(is, irep)=NaN;
         else
        [AUCs]=reconstruction_3papprox(round(x_all(:,2:end-1),3), ...
                    round(dt_x_all_v2,3), adjacency,num_nodes);
        all_AUCS(is, irep)=AUCs;
        end
    end
end

save(sprintf('simulations/summary_3p_round_sigma%d-%d_dt%0.0e_res%0.0e',...
    sigmamin, sigmamax, delta_t, res), ...
    'all_AUCS')
%% collect results for 3p alg, add noise
all_AUCS = zeros(Nnoise, Ntotal);

for is=1:Nnoise
    for irep=1:Ntotal
        load(sprintf(strcat(save_string, "_S%d_I%d.mat"), is, irep)); 
        x_all=x_all(:,2:end-1); 
        if sum(isnan(x_all(:)))>0
            all_AUCS(is, irep)=NaN;
         else
        [AUCs]=reconstruction_3papprox(x_all + sigma_add*randn(size(x_all)), ...
                    dt_x_all_v2 + sigma_add*randn(size(dt_x_all_v2)), ...
                    adjacency,num_nodes);
        all_AUCS(is, irep)=AUCs;
        end
    end
end

save(sprintf('simulations/summary_3p_addnoise_sigma%d-%d_dt%0.0e_res%0.0e',...
    sigmamin, sigmamax, delta_t, res), ...
    'all_AUCS')
%% nextstep alg, rounding

all_AUCS = zeros(Nnoise, Ntotal);

for is=1:Nnoise
    for irep=1:Ntotal
        load(sprintf(strcat(save_string, "_S%d_I%d.mat"), is, irep)); 
        if sum(isnan(x_all(:)))>0
            all_AUCS(is, irep)=NaN;
         else
        [AUCs]=reconstruction_nextstep(round(x_all,3), ...
                                    adjacency,num_nodes);
        all_AUCS(is, irep)=AUCs;
        end
    end
end

save(sprintf('simulations/summary_nextstep_round_sigma%d-%d_dt%0.0e_res%0.0e',...
    sigmamin, sigmamax, delta_t, res), ...
    'all_AUCS')

%% nextstep alg, add noise

all_AUCS = zeros(Nnoise, Ntotal);

for is=1:Nnoise
    for irep=1:Ntotal
        load(sprintf(strcat(save_string, "_S%d_I%d.mat"), is, irep));  
        if sum(isnan(x_all(:)))>0
            all_AUCS(is, irep)=NaN;
         else
        [AUCs]=reconstruction_nextstep(x_all + sigma_add*randn(size(x_all)), ...
                                    adjacency,num_nodes);
        all_AUCS(is, irep)=AUCs;
        end
    end
end
save(sprintf('simulations/summary_nextstep_addnoise_sigma%d-%d_dt%_dt%0.0e_res%0.0e',...
    sigmamin, sigmamax, delta_t, res), ...
    'all_AUCS')


%% plotting
pre_s=linspace(sigmamin,sigmamax,Nnoise);
pre_s= power(10,pre_s);  %log-spacing of noise


