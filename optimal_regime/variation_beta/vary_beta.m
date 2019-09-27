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



%lower res, same dt
 delta_t=0.001; %fine simulation tscale
 res=0.01;       %coarse simulations scale, used for reconstruction
beta
simulate=true;
 simulate_reconstruct_SNR(Nnoise, sigmamin, sigmamax,sigma_add,...
                        Ntotal,num_simulations,couplings,num_nodes,...
                          n_incoming, T, delta_t, res,initial, simulate)



function simulate_reconstruct_SNR(Nnoise, sigmamin, sigmamax,sigma_add,...
                       Ntotal,num_simulations,couplings,num_nodes,...
                         n_incoming, T, delta_t, res,initial, simulate)
%given parameters
%this function first simulated tseries witha  cubic decay
%reconstructs the connectivity matrices for both rounding and additional
%noise
%calculates SNR ration
disp(sigma_add);
save_string=sprintf('simulations/data_cubic_T%d_dt%0.0e_Nnodes%d_Nincoming%d_sigma%d-%d',...
    T,delta_t, num_nodes, n_incoming, sigmamin, sigmamax);

pre_s=linspace(sigmamin,sigmamax,Nnoise);
pre_s= power(10,pre_s);  %log-spacing of noise
%% simulate timeseries

if simulate
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


%% collect SNR, rounding
all_sd = zeros(Nnoise, Ntotal, num_nodes);

for is=1:Nnoise 
    var=pre_s(is)^2;
    for irep=1:Ntotal
        load(sprintf(strcat(save_string, "_S%d_I%d.mat"), is, irep)); 
        if sum(isnan(x_all(:)))>0
            all_sd(is, irep,:)=NaN;
         else
        x_all_noise = round(x_all(:,int16(end/2):end), 3);
        sds = (std(x_all_noise,0,2).^2 + mean(x_all_noise, 2))/var;
        all_sd(is, irep,:)=sds;
        end
    end
end

save(sprintf('simulations/summary_std_sigma_round_sigma%d-%d_dt%0.0e_res%0.0e',...
    sigmamin, sigmamax, delta_t, res), ...
    'all_sd')

%% collect SNR, add noise
all_sd = zeros(Nnoise, Ntotal, num_nodes);

for is=1:Nnoise 
    var=pre_s(is)^2;
    for irep=1:Ntotal
        load(sprintf(strcat(save_string, "_S%d_I%d.mat"), is, irep));  
        if sum(isnan(x_all(:)))>0
            all_sd(is, irep,:)=NaN;
         else
        x_all_noise = x_all(:,int16(end/2):end);
        x_all_noise = x_all_noise + sigma_add*randn(size(x_all_noise));
        sds = (std(x_all_noise,0,2).^2 + mean(x_all_noise, 2))/var;
        all_sd(is, irep,:)=sds;
        end
    end
end

save(sprintf('simulations/summary_std_sigma_addnoise_sigma%d-%d_dt%0.0e_res%0.0e',...
    sigmamin, sigmamax, delta_t, res), ...
    'all_sd')

%% collect SNR, full precision
all_sd = zeros(Nnoise, Ntotal, num_nodes);

for is=1:Nnoise 
    var=pre_s(is)^2;
    for irep=1:Ntotal
        load(sprintf(strcat(save_string, "_S%d_I%d.mat"), is, irep));  
        if sum(isnan(x_all(:)))>0
            all_sd(is, irep,:)=NaN;
         else
        x_all_noise = x_all(:,int16(end/2):end);        
        sds = (std(x_all_noise,0,2).^2 + mean(x_all_noise, 2))/var;
        all_sd(is, irep,:)=sds;
        end
    end
end

save(sprintf('simulations/summary_std_sigmafullprec_sigma%d-%d_dt%0.0e_res%0.0e',...
    sigmamin, sigmamax, delta_t, res), ...
    'all_sd')

end


