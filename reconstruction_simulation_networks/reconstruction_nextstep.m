function [AUCs]=reconstruction_nextstep(x_all, adjacency, num_nodes)
steady_state_all=zeros(num_nodes,1);
df_reconstr=deriv_steadystate_reconstr(x_all);
bin_adjacency = adjacency;
bin_adjacency(bin_adjacency>0)=1;
%bin_adjacency(logical(eye(size(bin_adjacency)))) = 1.;

bin_adjacency= reshape(bin_adjacency, num_nodes*num_nodes,1);
df_reconstr=reshape(df_reconstr, num_nodes*num_nodes,1);                                         
[~,~,~,AUCs] = perfcurve(bin_adjacency,...
                abs(df_reconstr), 1);



    function [df_reconstr]=deriv_steadystate_reconstr(x_all)
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

end