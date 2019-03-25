function [AUCs]=AUC_reconstr(x_all, dt_x_all, adjacency, num_nodes)
%calculate area under curve for binary reconstruction
steady_state_all=zeros(num_nodes,1);
df_reconstr=deriv_steadystate_reconstr(x_all, dt_x_all, steady_state_all);
bin_adjacency = adjacency;
bin_adjacency(bin_adjacency>0)=1;
AUCs=zeros(num_nodes,1);
for inode=1:num_nodes
    [fpr,tpr,~,AUC] = perfcurve(bin_adjacency(inode,1:end),...
                abs(df_reconstr(inode,1:end)), 1);
    AUCs(inode,1)=AUC ; 
end
end