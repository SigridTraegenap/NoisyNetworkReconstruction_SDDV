function [df_reconstr]=deriv_steadystate_reconstr(x_all, dt_x_all, steady_state_all)
%reconstruct possible connectivites (inouts) for single nodes
%use ML-estimates of weights contributing to target (dt_activity) and
%possible inputs (activity_mat)
%Moore-Penrose Pseudo-Inverse to calculate ML estimate
%ignore self-connectivity (decay terms), set to zero in both reconstruction
%and actual values
[num_nodes, num_steps]=size(x_all);
df_reconstr=nan(num_nodes, num_nodes);
for inode=1:num_nodes
    dt_activity=dt_x_all(inode,:);
    steady_state=repmat(steady_state_all(inode), 1,num_steps);
    activity_mat=x_all - ones(num_nodes,1)*steady_state;    
    %disp([isnan(activity_mat)]);
    act_PINV = pinv(activity_mat);     
    df_reconstr(inode,:)=mtimes(dt_activity,act_PINV);
end
df_reconstr(logical(eye(size(df_reconstr)))) = 0.;
end