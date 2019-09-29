function [df_reconstr]=Nextstep_deriv_steadystate_reconstr(x_all)
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