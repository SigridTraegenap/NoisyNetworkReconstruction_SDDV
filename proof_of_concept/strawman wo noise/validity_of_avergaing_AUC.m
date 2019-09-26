clear;

%load tseries
load(strcat('/Users/straegenap/Documents/Natwiss_Kolleg/NoisyNetworkReconstruction_SDDV/',...
'proof_of_concept/strawman wo noise/simulations/data_linear_T100_dt1e-03_',...
'Nnodes30_Nincoming5_sigma1e-01w_noise_I35.mat'))

%test AUC mean
[AUCs]=OLD_reconstruction_2p_approx(round(x_tau_all,3), ...
                    round(dt_x_all,3), adjacency,num_nodes);

disp([AUCs]);
disp(['mean AUC']);
disp(mean(AUCs));
%test AUC concatenate all 
%create df_reconstructed
[df_reconstr]=deriv_steadystate_reconstr(round(x_tau_all,3),...
                                             round(dt_x_all,3) );
bin_adjacency = adjacency;
bin_adjacency(bin_adjacency>0)=1; 
bin_adjacency= reshape(bin_adjacency, num_nodes*num_nodes,1);

df_reconstr=reshape(df_reconstr, num_nodes*num_nodes,1);
                                         
[~,~,~,AUC] = perfcurve(bin_adjacency,...
                abs(df_reconstr), 1);
 disp('concatenate AUC');
 disp(AUC);           