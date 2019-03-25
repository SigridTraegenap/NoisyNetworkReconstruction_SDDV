clear;
close all;

%simulate noisy time series in a sparse, coupled network

%add all higher order code
addpath("..\reconstruction_simulation_networks\")



%define parameter
num_nodes=30; %number_nodes
n_incoming=5; %number incoming nodes (sparsity)

delta_t=0.0001;  %fine simulation tscale
res=0.01;  %coarse reconstruction t-scale
initial="zero"; %specify initial conditions
T=100;  %how many timesteps are evaluated
couplings=[1];  %coupling strength 
num_simulations=1;  %how many simulations to perform 
%(only important when reconstructing w/o noise)

%where to save data
save_string=strcat("simulations1/data_linear_T100_dt-4_numnode"+num2str(num_nodes)+...
    "n_incoming"+num2str(n_incoming)+ "sigma-5+3");

pre_s=linspace(-5,3,100);
pre_s= power(10,pre_s);  %log-spacing of noise


for ia=1:size(couplings, 2)
    alpha=couplings(ia);
    for is=1:size(pre_s,2)
            sigma=pre_s(is); 
            disp([is, sigma]);
            [x_all, x_tau_all, dt_x_all, dt_x_all_v2, adjacency, betas]  = simulate_lineardecay_tseries(...
                sigma,num_simulations,alpha,num_nodes,...
                n_incoming, T, delta_t, res,initial);
            save(strcat(save_string+num2str(is)+...
                "alpha"+num2str(alpha)+".mat"))     
            
    end      
end



