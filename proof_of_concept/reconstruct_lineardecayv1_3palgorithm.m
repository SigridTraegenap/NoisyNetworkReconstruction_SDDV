clear;
close all;

%reconstruct connectivities with three point approx algorithm
%linear decay data
%intermediate approximation level

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

%show reconstructability with different time-lengths
times = linspace(0.1,1,10);
res_auc = zeros(size(times,2), size(pre_s,2));

for ia=1:size(couplings, 2)
    alpha=couplings(ia);
for is=1:size(pre_s,2)
        sigma=pre_s(is);  
        load(strcat(save_string+num2str(is)+...
                "alpha"+num2str(alpha)+".mat"));  
        x_all=x_all(:,2:end-1); 
        %change size to make it compatible with dt_x_all_v2   
        %loads adjacency,     
        disp(["noise level",log(sigma)/log(10),"iteration number", is])
        Ntot = size(x_all,2)-1;
        %disp([size(dt_x_all), size(x_all)]);
        for it=1:size(times, 2)            
            x_sample = x_all(:,1:int32(Ntot*times(it)));      
            %use other derivative approximation
            dt_x_sample = dt_x_all_v2(:,1:int32(Ntot*times(it)));
            if sum(isnan(dt_x_all(:)))>0
               res_auc(it, is)=NaN;
            else
                [AUCs]=reconstruction_3papprox(round(x_sample,3), ...
                    round(dt_x_sample,3), adjacency,num_nodes);
                res_auc(it, is)=mean(AUCs);
            end
        end        
        disp(["AUC for T=100",res_auc(end, is)]);
        drawnow('update');
end
end

save(strcat("simulations1/summmary_data_linear_T100_dt-4_numnode"+num2str(num_nodes)+...
    "n_incoming"+num2str(n_incoming)+ "sigma-5+3_3p"), "res_auc")


%preliminary plotting
figure()
imagesc(res_auc, [0.0 1])
c = colorbar;
c.Label.String = 'AUC';
ylabel("fraction of total time T=100")

Nstepy = int8(size(times,2)/6);
yticks(1:Nstepy:size(times,2))
yticklabels(times(1:Nstepy:size(times,2)))

xlabel("noise")
Nstepx = int8(size(pre_s,2)/5);
xticks(1:Nstepx:size(pre_s,2))
xticklabels(pre_s(1:Nstepx:size(pre_s,2)))
key=strcat("linear_3p_paramscan_sigma-5+3_alpha"+num2str(alpha));
print(key, "-dpdf")