function [betas,adjacency]=setup_NW(num_nodes, ni, coupling, beta)
%set up adjacency matrix of network
% sparse connectivity (ni incoming nodes)
% coupling is a fixed parameter, thus all coupled nodes are coupled
% with a similar strength
adjacency=zeros(num_nodes, num_nodes);
fiedler_eig=-1;
while fiedler_eig<0
    for inode=1:num_nodes
        f=randperm(num_nodes);
        f=f(1:ni);
        adjacency(inode,f) = coupling;
    end
    adjacency(logical(eye(size(adjacency)))) = 0.;
bin_adjacency= adjacency;
bin_adjacency(bin_adjacency>0)=1;
e = abs(eig(bin_adjacency));
Xs = sort(e);
fiedler_eig = Xs(2,1);
if fiedler_eig<0
    disp("not connected");
end
end
betas = ones(num_nodes,1)*beta;
end