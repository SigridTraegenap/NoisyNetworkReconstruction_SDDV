function [AUCs]=reconstruction_3papprox(x_all, dt_x_all_v2, adjacency, num_nodes)
steady_state_all=zeros(num_nodes,1);
%x_all = x_all(:,2:end-1);
%input to the functions should have same length 
%need to cut input x_all in outer function
df_reconstr=deriv_steadystate_reconstr(x_all, dt_x_all_v2, steady_state_all);
bin_adjacency = adjacency;
bin_adjacency(bin_adjacency>0)=1;
bin_adjacency(logical(eye(size(bin_adjacency)))) = 1.;
AUCs=zeros(num_nodes,1);
for inode=1:num_nodes
    [~,~,~,AUC] = perfcurve(bin_adjacency(inode,1:end),...
                abs(df_reconstr(inode,1:end)), 1);
    AUCs(inode,1)=AUC ; 
end


function [df_reconstr]=deriv_steadystate_reconstr(x_tau, dt_x_all, steady_state_all)
[num_nodes, num_steps]=size(x_tau);
df_reconstr=nan(num_nodes, num_nodes);
for inode2=1:num_nodes
    dt_activity=dt_x_all(inode2,:);
    steady_state=repmat(steady_state_all(inode2), 1,num_steps);
    activity_mat=x_tau - ones(num_nodes,1)*steady_state;    
    act_PINV = pinv(activity_mat);     
    df_reconstr(inode2,:)=mtimes(dt_activity,act_PINV);
end
df_reconstr(logical(eye(size(df_reconstr)))) = 0.;
end

end