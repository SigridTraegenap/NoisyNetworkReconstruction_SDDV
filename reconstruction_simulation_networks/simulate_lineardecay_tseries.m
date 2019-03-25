function [x_all, x_tau_all, dt_x_all, dt_x_all_v2, adjacency, betas] = simulate_lineardecay_tseries(...
    sigma,num_simulations,coupling, num_nodes, n_incoming, T, delta_t,...
    res,initial_cond)
Ncoarse = T/res -1;
[betas,adjacency]=setup_NW(num_nodes, n_incoming, coupling, 1);
x_all=nan(num_nodes,num_simulations*(Ncoarse-1));
dt_x_all=nan(num_nodes,num_simulations*(Ncoarse-1));
x_tau_all=nan(num_nodes,num_simulations*(Ncoarse-1));
dt_x_all_v2=nan(num_nodes,num_simulations*(Ncoarse-2));
for isim=1:num_simulations
    [~,activitiyC]=simulate_networknoise(betas,adjacency,T,delta_t, res, sigma,initial_cond);     
    [~, x_tau, derivative]=time_derivative_approx(activitiyC, T, res);  
    derivative2=time_derivative_IMapprox(activitiyC, res);
    x_all(:,(isim-1)*(Ncoarse)+1:isim*(Ncoarse)+1)=activitiyC;
    x_tau_all(:,(isim-1)*(Ncoarse-1)+1:isim*(Ncoarse-1)+1)=x_tau;
    dt_x_all(:,(isim-1)*(Ncoarse-1)+1:isim*(Ncoarse-1)+1)=derivative;
    dt_x_all_v2(:,(isim-1)*(Ncoarse-2)+1:isim*(Ncoarse-2)+1)=derivative2;
end

    function [activityF, activityC]=simulate_networknoise(betas,adjacency,T,...
        delta_t, res, sigma,initial_cond )
    num_nodes=size(betas,1);
    num_steps=T/delta_t;
    sqrtdt =sqrt(delta_t);
    ind_sampling = 1:int32(res/delta_t):int32(num_steps);
    activityF=zeros(num_nodes, num_steps+1);

    if contains(initial_cond, "random")
        activityF(:,1)=rand(num_nodes,1);
    end 

    for istep=1:num_steps 
        S= rand(num_nodes,1);
        S(S>=0.5)=1;
        S(S<0.5)=-1;
        W=sqrtdt.*randn(num_nodes,1);
        K1 = delta_t*update_nw(activityF(:,istep), adjacency, betas) +sigma.*(W-sqrtdt*S);
        K12 = activityF(:,istep) + K1; 
        K2 = delta_t*update_nw(K12, adjacency, betas) +sigma.*(W+sqrtdt*S);
        activityF(:,istep+1) = activityF(:,istep)+0.5*(K1+K2);
    end
    activityC =activityF(:,ind_sampling);

        function [dy]=update_nw(act, adjacency,betas)
        outgoing_adj=sum(adjacency,2);
        dy= -betas.*power(act,1) +...
                adjacency*act -...
                outgoing_adj.*act;
        end
    end

end

