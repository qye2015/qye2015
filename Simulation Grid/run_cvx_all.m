% run cvx for different approaches:
%
% 1. BS Arch
% 2. Mix Arch
% 3. Same-size Arch
% 4. Orthogonal, BS arch
% 5. Orthogonal, cellular case
% 6. max-SINR
% 7. LB for cellular case
% 8. eICIC
% 
% Network deployment: 
%   1. Small cells are uniformly dropped
%   2. Macro BS is located at the center of the network
%   3. Users are uniformly distributed with density nu
%   4. All users are static for now

% Last updated: 1/29/14 3:30pm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; clf; close all;
% ------------------------------------------
step_rho                            =   1; % 0.75 0.5 0.25];  
% how many nearby BSs may become a cluster
N_nearBS                            =   8;
N_nearBS_orth                       =   7;
% ------------------------------------------
load(strcat('data/',int2str(1),'-nphot-',int2str(3),'-nuhot-',int2str(6),'-grid.mat'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% run data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
loop                =   0;
save('data/loop.mat','loop');
    
totloop             =   1;
while loop < totloop
        
        load('data/loop.mat');
        loop                        =   loop+1;
        if mod(loop, 4) ==0
            clc;
        end
        fprintf('start loop=%d\n',loop);
        
        load(strcat('data/',int2str(loop),'-nphot-',int2str(denseP_hot/denseP_reg),'-nuhot-',int2str(denseu_hot/denseu_reg),'-grid.mat'));
        load(strcat('data/','rate-',int2str(loop),'-rho-0',int2str(step_rho*100),'-nphot-',int2str(denseP_hot/denseP_reg),'-nuhot-',int2str(denseu_hot/denseu_reg),'-grid.mat'));


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% cvx for mix arch, where clusters with different sizes can work together
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % number of clusters per user
%         nCluster_perU           =   size(act_u_cluster,2);
%        
%         
%         % =================================================
%         % the matrix to get j\in C, |C|=s
%         actBS_arch              =   zeros(nB*narch, nCluster);
%         for i = 1 : narch
%             actBS_arch(nB*(i-1)+1:nB*i, :) = actBSv'.*repmat(archv(i,:),nB, 1);
%         end
%         
%         % the matrix indicating j\in C |C|=s for each user
%         % col:(j1s1, j2s1, ..., j1s2, j2s2, ...)
%         % row: c1u1, c2u1, ..., c1u2, c2u2, ...
%         actBS_arch_u            =   zeros(nB * narch, nCluster_perU * nu);
%         for i = 1 : nCluster_perU
%             actBS_arch_u(:, i:nCluster_perU:nu*nCluster_perU) = actBS_arch(:, act_u_cluster(:,i));
%         end
%         
%         % =================================================
%         % the matrix to get |C|=s
%         actcluster_arch_u     	=   zeros(narch*nu, nCluster_perU);
%         for i = 1 : narch
%             actcluster_arch_u(nu*(i-1)+1: nu*i, :) = reshape(archv(i, act_u_cluster),nu,nCluster_perU);
%         end
        % =================================================
        % MIMO parameters                       
        % user factor
        rho                         =   ones(nB, 1);
        for i = 2 : nB
            rho(i)                  =   step_rho * i;
        end
        % =================================================
        % number of users that can be served
        SU                      =   ones(nM + nP, 1);
        SU(1 : nM)              =   S_macro;
        SU(nM+1 : nM+nP)     	=   S_pico;
        SU_js                   =   zeros(nB*narch, 1);
        for i = 1 : narch
            SU_js(nB*(i-1)+1:nB*i, :) = round(max(SU, SU*rho(i)));
        end
        
%         % matrix to indicator 
%         % row: j1s1(u1, u2, u3, ...), j2s1(u1,u2,u3,...), ..., j1s2(u1,u2,u3...), j2s2(u1,u2,u3)
%         % col: c1, c2, c3, ...
%         act_arch_u_cluster      =   reshape((reshape(actBS_arch_u,[],nu))',[], nCluster_perU);
%         % =================================================
% %         for i = 1:nu
% %             c_zf_rho_2(:, i) = c_zf_rho(act_u_cluster(i,:),i);
% %         end
%         % =================================================
%         % start cvx 
%         start_time = tic;
%         cvx_expert true
%         cvx_clear
%         cvx_begin
%             variable y(nCluster_perU, nu);
%             variable mu(nB, narch);
%             maximize( sum(log(sum(c_zf_rho.*y, 1))) );
%             %maximize( geo_mean(sum(c_zf_rho.*y, 1)) );
%         subject to
%             actBS_arch_u*reshape(y, nCluster_perU*nu, 1)./SU_js <= reshape(mu,nB*narch,1);
%             sum(act_arch_u_cluster.*repmat(y',nB*narch,1),2) <= kron(reshape(mu,nB*narch,1), ones(nu,1));
%             sum(y, 1) <= 1;
%             sum(mu, 2) <= 1;
%             y >= 0;
%             mu >= 0;
%         cvx_end;
%         cvx_time_mixA           =   toc(start_time);
%         cvx_status_mixA         =   cvx_status;
% 
%         % =================================================
%         % data post-process
%         y2                      =   y;
%         y2(y<10^(-5))           =   0;
%         y2(y>1-10^(-5))         =   1;
%         y_cvx_mixA              =   y2;
%         mu2                     =   mu;
%         mu2(mu<10^(-5))         =   0;
%         mu2(mu>1-10^(-5))       =   1;
%         mu_cvx_mixA             =   mu2;
% 
%         % =================================================
%         % find fractional users -- users connecting to multiple clusters in an arch.
%         ytp                     =   y_cvx_mixA;
%         ytp(y_cvx_mixA>0)       =   1;
%         % how many clusters with each size the user connecting to: nA*nu
%         ytp_perS                =   (reshape(sum(actcluster_arch_u .* repmat(ytp', narch, 1),2), nu, narch))';
%         % -------------------------------------------
%         % find the total number of fractional users
%         frac_ue_id              =   find(max(ytp_perS, [], 2) > 1);
%         frac_ue_mixA            =   length(frac_ue_id);
%         % find the number of fractional users per s (i.e., clusters with size s)
%         frac_ue_indic_perS      =   ytp_perS > 1;
%         frac_ue_perS_mixA       =   sum(frac_ue_indic_perS, 2);
%         % =================================================
%         % matrix to indicate macro BS
%         macro_index             =   find(sum(actBSv(:,1:nM),2)>0 & sum(actBSv,2) ==1);
%         actcluster_macro_u      =   zeros(nCluster_perU, nu);
%         for i = 1 : length(macro_index)
%             actcluster_macro_u   = 	min(1, actcluster_macro_u + (act_u_cluster == macro_index(i))');
%         end
%         % find users only connecting to macro BSs
%         if ~isempty(macro_index)
%             ytp_macro           =   ytp .* actcluster_macro_u;
%             macro_uetp          =   find(max(ytp_macro,[],1)>0 & max(ytp_perS, [], 1) <= 1);
%             macro_ue_mixA       =   length(macro_uetp);
% 
%             % users that is served by macro BSs in any instance
%             macro_ue_only_mixA	=   length(find(max(ytp_macro,[],1)>0 & sum(ytp,1)<=1));
%         end
%         % =================================================
%         % matrix to indicate pico BS
%         pico_index              =   find(sum(actBSv(:,1+nM:nB),2)>0 & sum(actBSv,2) ==1);
%         actcluster_pico_u       =   zeros(nCluster_perU, nu);
%         for i = 1 : length(pico_index)
%             actcluster_pico_u   = 	min(1, actcluster_pico_u + (act_u_cluster == pico_index(i))');
%         end
%         % find users only connecting to small BSs
%         if ~isempty(pico_index)
%             ytp_pico            =   ytp .* actcluster_pico_u;
%             pico_uetp           =   find(max(ytp_pico,[],1)>0 & max(ytp_perS, [], 1) <= 1);
%             pico_ue_mixA        =   length(pico_uetp);
% 
%             % users that is served by macro BSs in any instance
%             pico_ue_only_mixA	=   length(find(max(ytp_pico,[],1)>0 & sum(ytp,1)<=1));
%         end
%         % =================================================
%         % matrix to indicate clusters
%         cluster_index          	=   find(sum(actBSv,2) > 1);
%         actcluster_cluster_u 	=   zeros(nCluster_perU, nu);
%         for i = 1 : length(cluster_index)
%             actcluster_cluster_u=   min(1, actcluster_cluster_u + (act_u_cluster == cluster_index(i))');
%         end
%         % find users only connecting to clusters with multiple BSs
%         if ~isempty(cluster_index)
%             ytp_cluster       	=   ytp .* actcluster_cluster_u;
%             cluster_uetp     	=   find(max(ytp_cluster,[],1)>0 & max(ytp_perS, [], 1) <= 1);
%             cluster_ue_mixA    	=   length(cluster_uetp);
% 
%             % users that is served by macro BSs in any instance
%             cluster_ue_only_mixA=   length(find(max(ytp_cluster,[],1)>0 & sum(ytp,1)<=1));
%         end
% 
%         % =================================================
%         % find how long each BS working in clusters
%         BS_time_mixA            =   mu_cvx_mixA;
% 
%         % =================================================
%         % utility function
%         rate_cvx_mixA           =   sum(c_zf_rho.*y_cvx_mixA, 1);
%         utility_cvx_mixA        =   sum(log(rate_cvx_mixA));
% 
        [y_cvx_mixA, mu_cvx_mixA, cvx_time_mixA, cvx_status_mixA, frac_ue_mixA, frac_ue_perS_mixA, macro_ue_mixA,  macro_ue_only_mixA, pico_ue_mixA, pico_ue_only_mixA, cluster_ue_mixA, cluster_ue_only_mixA, BS_time_mixA, rate_cvx_mixA, utility_cvx_mixA] = cvx_mixBSarch(nu, nM, nP, nB, narch, nCluster, SU_js, act_u_cluster, actBSv, archv, c_zf_rho);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% cvx for predefined arch
        % clusters with size 1 are allowed in all architectures
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % add 
%         nCluster_perU_predef    =   nCluster_perU + nB * (narch - 1);
%         
%         % =================================================
%         % Each predefined arch includes clusters with same size and individual BSs
%         % rate matrix in the predef arch
%         c_zf_predef             =   [c_zf_rho; repmat(c_zf_BS, 3, 1)];
%         actBS_arch_u_predef     =   [zeros(nB, nB*(narch-1)*nu); ones(nB, nB*nu), zeros(nB, nB*(narch-2)*nu); zeros(nB,nB*nu), ones(nB,nB*nu), zeros(nB,nB*nu); zeros(nB, nB*(narch-2)*nu), ones(nB, nB*nu)];
%         actcluster_arch_u_predef=   [zeros(1,nB*(narch-1));ones(1,nB), zeros(1,nB*(narch-2));zeros(1,nB),ones(1,nB),zeros(1,nB);zeros(1,nB*(narch-2)),ones(1,nB)];
%         % =================================================
%         % start cvx 
%         start_time = tic;
%         cvx_expert true
%         cvx_clear
%         cvx_begin
%             variable y(nCluster_perU_predef, nu);
%             variable lambda(narch, 1);
%             maximize( sum(log(sum(c_zf_predef.*y, 1))) );
%             %maximize( geo_mean(sum(c_zf_rho.*y, 1)) );
%         subject to
%             actBS_arch_u*reshape(y(1:nCluster_perU,:), nCluster_perU*nu, 1)./SU_js +  actBS_arch_u_predef*reshape(y(1+nCluster_perU:nCluster_perU_predef,:), nB*(narch-1)*nu, 1)./SU_js <= kron(lambda,ones(nB,1));
%             sum(actcluster_arch_u.*repmat((y(1:nCluster_perU,:))',narch,1),2) + reshape((actcluster_arch_u_predef*y(nCluster_perU+1:nCluster_perU_predef,:))',[],1) <= kron(lambda, ones(nu,1));
%             %sum(y, 1) <= 1;
%             sum(lambda) <= 1;
%             y >= 0;
%             lambda >= 0;
%         cvx_end;
%         cvx_time_predef         =   toc(start_time);
%         cvx_status_predef       =   cvx_status;
%         
%         % =================================================
%         % data post-process
%         y2                      =   y;
%         y2(y<10^(-5))           =   0;
%         y2(y>1-10^(-5))         =   1;
%         y_cvx_predef           	=   y2;
%         lambda2              	=   lambda;
%         lambda2(lambda<10^(-5))	=   0;
%         lambda2(lambda>1-10^(-5))=   1;
%         lambda_cvx_predef    	=   lambda2;
% 
%         % =================================================
%         % find fractional users -- users connecting to multiple clusters in an arch.
%         ytp                     =   y_cvx_predef;
%         ytp(y_cvx_predef>0) 	=   1;
%         % how many clusters (no individual BS) with each size the user connecting to : nA * nu
%         ytp_perS_predef_c   	=   (reshape(sum(actcluster_arch_u .* repmat(ytp(1:nCluster_perU,:)', narch, 1),2), nu, narch))';
%         % how many clusters (only individual BS) with each size the user connecting to : nA * nu
%         ytp_perS_predef_BS      =   actcluster_arch_u_predef*ytp(nCluster_perU+1:nCluster_perU_predef,:);
%         % -------------------------------------------
%         % find the total number of fractional users
%         frac_ue_id_predef    	=   find(max(max(ytp_perS_predef_c, [], 2),max(ytp_perS_predef_BS, [], 2)) > 1);
%         frac_ue_predef       	=   length(frac_ue_id_predef);
%         % find the number of fractional users per s (i.e., clusters with size s)
%         frac_ue_indic_perS_predef=   (ytp_perS_predef_c > 1 | ytp_perS_predef_BS > 1);
%         frac_ue_perS_predef    	=   sum(frac_ue_indic_perS_predef, 2);
%         
%         % =================================================
%         % matrix to indicate macro BS
%         macro_index_pf          =   [];
%         for i = 1 : (narch-1)
%             macro_index_pf      =   [macro_index_pf, nCluster_perU + nB*(i-1)+[1:nM]];
%         end
%         % find users only connecting to macro BSs
%         if ~isempty(macro_index)
%             ytp_macro_pf_c      =   ytp(1:nCluster_perU,:) .* actcluster_macro_u;
%             macro_uetp          =   find((max(ytp_macro_pf_c,[],1)>0 | sum(ytp(macro_index_pf,:),1)>0) & max(ytp_perS_predef_c, [], 1) <= 1 & max(ytp_perS_predef_BS, [], 1) <= 1);
%             macro_ue_predef     =   length(macro_uetp);
% 
%             % users that is served by macro BSs in any instance
%             macro_ue_only_predef=   length(find((max(ytp_macro_pf_c,[],1)>0 | sum(ytp(macro_index_pf,:),1)>0) & (sum(ytp(1:nCluster_perU,:),1)<=1 | sum(ytp(nCluster_perU+1:nCluster_perU_predef,:),1)<=1)));
%         end
%         
%         % =================================================
%         % matrix to indicate pico BS
%         pico_index_pf           =   [];
%         for i = 1 : (narch-1)
%             pico_index_pf       =   [pico_index_pf, nCluster_perU + nB*(i-1)+[1+nM:nB]];
%         end
%         % find users only connecting to macro BSs
%         if ~isempty(pico_index)
%             ytp_pico_pf_c       =   ytp(1:nCluster_perU,:) .* actcluster_pico_u;
%             pico_uetp           =   find((max(ytp_pico_pf_c,[],1)>0 | sum(ytp(pico_index_pf,:),1)>0) & max(ytp_perS_predef_c, [], 1) <= 1 & max(ytp_perS_predef_BS, [], 1) <= 1);
%             pico_ue_predef      =   length(pico_uetp);
% 
%             % users that is served by macro BSs in any instance
%             pico_ue_only_predef =   length(find((max(ytp_pico_pf_c,[],1)>0 | sum(ytp(pico_index_pf,:),1)>0) & (sum(ytp(1:nCluster_perU,:),1)<=1 | sum(ytp(nCluster_perU+1:nCluster_perU_predef,:),1)<=1)));
%         end  
%         
%         % =================================================
%         % matrix to indicate cluster
%         % find users only connecting to macro BSs
%         if ~isempty(cluster_index)
%             ytp_cluster_pf_c  	=   ytp(1:nCluster_perU,:) .* actcluster_cluster_u;
%             Cluster_uetp      	=   find(max(ytp_Cluster_pf_c,[],1)>0  & max(ytp_perS_predef_c, [], 1) <= 1 & max(ytp_perS_predef_BS, [], 1) <= 1);
%             Cluster_ue_predef 	=   length(Cluster_uetp);
% 
%             % users that is served by macro BSs in any instance
%             Cluster_ue_only_predef=	length(find(max(ytp_cluster_pf_c,[],1)>0 & sum(ytp(1:nCluster_perU,:),1)<=1 ));
%         end
%         
%         % =================================================
%         % find how long each BS working in clusters
%         BS_time_predef          =   reshape(actBS_arch_u*reshape(y_cvx_predef(1:nCluster_perU, :), nCluster_perU*nu, 1)./SU_js, nB, narch);
%         for i = 1 : nB
%             indextp             =   nCluster_perU + [i:nB:nB*(narch-1)];
%             BS_time_predef(i,1) =   BS_time_predef(i, 1) + sum(sum(y_cvx_predef(indextp,:)))/SU_js(i);
%         end
%         
%         % =================================================
%         % utility function
%         rate_cvx_predef         =   sum(c_zf_predef.*y_cvx_predef, 1);
%         utility_cvx_predef  	=   sum(log(rate_cvx_predef));

        [y_cvx_predef, lambda_cvx_predef, cvx_time_predef, cvx_status_predef, frac_ue_predef, frac_ue_perS_predef, macro_ue_predef,  macro_ue_only_predef, pico_ue_predef, pico_ue_only_predef, Cluster_ue_predef, Cluster_ue_only_predef, BS_time_predef, rate_cvx_predef, utility_cvx_predef] = cvx_MixArch(nu, nM, nP, nB, narch, nCluster, SU_js, act_u_cluster, actBSv, archv, c_zf_rho, c_zf_BS);

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% cvx for time sharing of archs, where each arch has same-size clusters
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % =================================================
%         % Each predefined arch includes clusters with same size 
%         % start cvx 
%         start_time = tic;
%         cvx_expert true
%         cvx_clear
%         cvx_begin
%             variable y(nCluster_perU, nu);
%             variable lambda(narch, 1);
%             maximize( sum(log(sum(c_zf_rho.*y, 1))) );
%             %maximize( geo_mean(sum(c_zf_rho.*y, 1)) );
%         subject to
%             actBS_arch_u*reshape(y, nCluster_perU*nu, 1)./SU_js <= kron(lambda,ones(nB,1));
%             sum(actcluster_arch_u.*repmat(y',narch,1),2) <= kron(lambda, ones(nu,1));
%             %sum(y, 1) <= 1;
%             sum(lambda) <= 1;
%             y >= 0;
%             lambda >= 0;
%         cvx_end;
%         cvx_time_sizearch 	=   toc(start_time);
%         cvx_status_sizearch =   cvx_status;
%         % =================================================
%         % data post-process
%         y2                      =   y;
%         y2(y<10^(-5))           =   0;
%         y2(y>1-10^(-5))         =   1;
%         y_cvx_sizearch       	=   y2;
%         lambda2              	=   lambda;
%         lambda2(lambda<10^(-5))	=   0;
%         lambda2(lambda>1-10^(-5))=   1;
%         lambda_cvx_sizearch   	=   lambda2;
%         
%         % =================================================
%         % find fractional users -- users connecting to multiple clusters in an arch.
%         ytp                     =   y_cvx_sizearch;
%         ytp(y_cvx_sizearch>0) 	=   1;
%         % how many clusters (no individual BS) with each size the user connecting to : nA * nu
%         ytp_perS_sizearch   	=   (reshape(sum(actcluster_arch_u .* repmat(ytp, narch, 1),2), nu, narch))';
%         % -------------------------------------------
%         % find the total number of fractional users
%         frac_ue_id_sizearch    	=   find(max(ytp_perS_sizearch_c, [], 2) > 1);
%         frac_ue_sizearch       	=   length(frac_ue_id_sizearch);
%         % find the number of fractional users per s (i.e., clusters with
%         % size s)
%         frac_ue_indic_perS_sizearch=   ytp_perS_sizearch > 1;
%         frac_ue_perS_sizearch  	=   sum(frac_ue_indic_perS_sizearch, 2);
%         
%         % =================================================
%         % find users only connecting to macro BSs
%         if ~isempty(macro_index)
%             ytp_macro_sizearch 	=   ytp .* actcluster_macro_u;
%             macro_uetp          =   find(max(ytp_macro_sizearch,[],1)>0  & max(ytp_perS_sizearch, [], 1) <= 1);
%             macro_ue_sizearch  	=   length(macro_uetp);
%             % users that is served by macro BSs in any instance
%             macro_ue_only_sizearch=   length(find(max(ytp_macro_sizearch,[],1)>0 & sum(ytp,1)<=1 ));
%         end
%         
%         % =================================================
%         % find users only connecting to macro BSs
%         if ~isempty(pico_index)
%             ytp_pico_sizearch  	=   ytp .* actcluster_pico_u;
%             pico_uetp           =   find(max(ytp_pico_sizearch,[],1)>0 & max(ytp_perS_sizearch, [], 1) <= 1);
%             pico_ue_sizearch   	=   length(pico_uetp);
% 
%             % users that is served by macro BSs in any instance
%             pico_ue_only_sizearch =   length(find(max(ytp_pico_sizearch,[],1)>0 & sum(ytp,1)<=1 ));
%         end  
%         
%         % =================================================
%         % find users only connecting to macro BSs
%         if ~isempty(cluster_index)
%             ytp_cluster_sizearch=   ytp .* actcluster_cluster_u;
%             Cluster_uetp      	=   find(max(ytp_Cluster_sizearch,[],1)>0  & max(ytp_perS_sizearch, [], 1) <= 1);
%             Cluster_ue_sizearch =   length(Cluster_uetp);
%             % users that is served by macro BSs in any instance
%             Cluster_ue_only_sizearch=	length(find(max(ytp_cluster_sizearch,[],1)>0 & sum(ytp,1)<=1 ));
%         end
%         
%         % =================================================
%         % find how long each BS working in clusters
%         BS_time_sizearch     	=   reshape(actBS_arch_u*reshape(y_cvx_sizearch, nCluster_perU*nu, 1)./SU_js, nB, narch);
% 
%         % =================================================
%         % utility function
%         rate_cvx_sizearch       =   sum(c_zf_rho.*y_cvx_sizearch, 1);
%         utility_cvx_sizearch	=   sum(log(rate_cvx_sizearch));

        [y_cvx_sizearch, lambda_cvx_sizearch, cvx_time_sizearch, cvx_status_sizearch, frac_ue_sizearch, frac_ue_perS_sizearch, macro_ue_sizearch,  macro_ue_only_sizearch, pico_ue_sizearch, pico_ue_only_sizearch, Cluster_ue_sizearch, Cluster_ue_only_sizearch, BS_time_sizearch, rate_cvx_sizearch, utility_cvx_sizearch] = cvx_SamesizeArch(nu, nM, nP, nB, narch, nCluster, SU_js, act_u_cluster, actBSv, archv, c_zf_rho);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% orthogonal, cvx for BS Arch
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                     
        % =================================================
        % number of users that can be served in orthogonal case
        SU_js_orth                  =   zeros(nB*narch_orth, 1);
        SU_js_orth(1:nB)            =   SU;
        for i = 2 : narch_orth
            SU_js_orth(nB*(i-1)+1:nB*i, :) = round(max(SU, SU*rho(i-1)));
        end
        actBSv_orth(1:nM,1:nM)      =   eye(nM);
        [y_cvx_mixA_orth, mu_cvx_mixA_orth, cvx_time_mixA_orth, cvx_status_mixA_orth, frac_ue_mixA_orth, frac_ue_perS_mixA_orth, macro_ue_mixA_orth,  macro_ue_only_mixA_orth, pico_ue_mixA_orth, pico_ue_only_mixA_orth, cluster_ue_mixA_orth, cluster_ue_only_mixA_orth, BS_time_mixA_orth, rate_cvx_mixA_orth, utility_cvx_mixA_orth] = cvx_mixBSarch(nu, nM, nP, nB, narch_orth, nCluster_orth, SU_js_orth, act_u_cluster_orth, actBSv_orth, archv_orth, c_zf_orth);

        
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% orthogonal, cvx for time sharing of archs, where each arch has same-size clusters
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [y_cvx_sizearch_orth, lambda_cvx_sizearch_orth, cvx_time_sizearch_orth, cvx_status_sizearch_orth, frac_ue_sizearch_orth, frac_ue_perS_sizearch_orth, macro_ue_sizearch_orth,  macro_ue_only_sizearch_orth, pico_ue_sizearch_orth, pico_ue_only_sizearch_orth, Cluster_ue_sizearch_orth, Cluster_ue_only_sizearch_orth, BS_time_sizearch_orth, rate_cvx_sizearch_orth, utility_cvx_sizearch_orth] = cvx_SamesizeArch(nu, nM, nP, nB, narch_orth, nCluster_orth, SU_js_orth, act_u_cluster_orth, actBSv_orth, archv_orth, c_zf_orth);

        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% max-sinr (no cluster)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [x_max, act_frac_max, rate_max, utility_max] =    max_bias_SINR(2.^c_zf_BS-1, bias, c_zf_BS, SU);
        % data post-process
        x2                      =   x_max;
        x2(x_max<10^(-5))       =   0;
        x2(x_max>1-10^(-5))     =   1;
        % find users only connecting to macro BSs
        %xtp                    =   zeros(size(x2));
        %xtp(1:nM,:)            =   x2(1:nM,:);
        %macro_ue_max           =   find(sum(xtp, 1) > 0);
        % find users only connecting to pico BSs
        xtp                     =   zeros(size(x2));
        xtp(nM+1:nM+nP,:)       =   x2(nM+1:nM+nP,:);
        pico_ue_max             =   length(find(sum(xtp, 1) > 0));



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% begin cvx (no cluster)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % =================================================
        start_time = tic;
        cvx_expert true
        cvx_clear
        cvx_begin
            variable y(nB, nu);
            maximize( sum(log(sum(c_zf_BS.*y, 1))) );
        subject to
            sum(y, 2) <= SU;
            sum(y, 1) <= 1;
            y >= 0;
        cvx_end;
        cvx_elapsed_time_BS  =   toc(start_time);

        % =================================================
        % data post-process
        y2                      =   y;
        y2(find(y<10^(-5)))     =   0;
        y2(find(y>1-10^(-5)))   =   1;
        x_cvx_BS                =   y2;
        % find fractional users
        ytp                     =   y2;
        ytp(y2>0)               =   1;
        frac_ue_BS_index        =   find(sum(ytp, 1) > 1);
        frac_ue_BS              =   length(frac_ue_BS_index);
        % utility function
        rate_cvx_BS             =   sum(c_zf_BS.*y2, 1);
        utility_cvx_BS          =   sum(log(rate_cvx_BS));

        % find users only connecting to macro BSs
        % xtp                 =   zeros(size(y2));
        %xtp(1:nM,:)         =   y2(1:nM,:);
        %xtp(:,frac_ue_BS)   =   0;
        %macro_ue_BS         =   find(sum(xtp, 1) > 0);
        % find users only connecting to pico BSs
        xtp                     =   zeros(size(y2));
        xtp(nM+1:nM+nP,:)       =   y2(nM+1:nM+nP,:);
        if frac_ue_BS > 0
            xtp(:,frac_ue_BS_index)	=   0;
        end
        pico_ue_BS              =   length(find(sum(xtp, 1) > 0));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% eICIC
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [y_cvx_ic, lambda_cvx_ic, cvx_time_ic, cvx_status_ic, frac_ue_ic, frac_ue_perS_ic, macro_ue_ic,  macro_ue_only_ic, pico_ue_ic, pico_ue_only_ic, Cluster_ue_ic, Cluster_ue_only_ic, BS_time_ic, rate_cvx_ic, utility_cvx_ic] = cvx_eICIC(nu, nM, nP, nB, nCluster, narch, SU_js, SU, act_u_cluster, actBSv, archv, c_zf_rho, c_zf_off);


        %utility_v               =   [utility_cvx_mixA, utility_cvx_mixA_ub, utility_cvx_predef, utility_cvx_arch2, utility_cvx_BS, utility_max_m];
        save(strcat('data/','cvx-',int2str(loop),'-rho-0',int2str(step_rho*100),'-nphot-',int2str(denseP_hot/denseP_reg),'-nuhot-',int2str(denseu_hot/denseu_reg),'-grid.mat'), 'y_cvx_mixA','mu_cvx_mixA','cvx_time_mixA','cvx_status_mixA','frac_ue_mixA','frac_ue_perS_mixA','macro_ue_mixA','macro_ue_only_mixA','pico_ue_mixA','pico_ue_only_mixA','cluster_ue_mixA','cluster_ue_only_mixA','BS_time_mixA','rate_cvx_mixA','utility_cvx_mixA','y_cvx_predef','lambda_cvx_predef','cvx_time_predef','cvx_status_predef','frac_ue_predef','frac_ue_perS_predef','macro_ue_predef','macro_ue_only_predef','pico_ue_predef','pico_ue_only_predef','Cluster_ue_predef','Cluster_ue_only_predef','BS_time_predef','rate_cvx_predef','utility_cvx_predef','y_cvx_sizearch','lambda_cvx_sizearch','cvx_time_sizearch','cvx_status_sizearch','frac_ue_sizearch','frac_ue_perS_sizearch','macro_ue_sizearch','macro_ue_only_sizearch','pico_ue_sizearch','pico_ue_only_sizearch','Cluster_ue_sizearch','Cluster_ue_only_sizearch','BS_time_sizearch','rate_cvx_sizearch','utility_cvx_sizearch','y_cvx_mixA_orth','mu_cvx_mixA_orth','cvx_time_mixA_orth','cvx_status_mixA_orth','frac_ue_mixA_orth','frac_ue_perS_mixA_orth','macro_ue_mixA_orth','macro_ue_only_mixA_orth','pico_ue_mixA_orth','pico_ue_only_mixA_orth','cluster_ue_mixA_orth','cluster_ue_only_mixA_orth','BS_time_mixA_orth','rate_cvx_mixA_orth','utility_cvx_mixA_orth','y_cvx_sizearch_orth','lambda_cvx_sizearch_orth','cvx_time_sizearch_orth','cvx_status_sizearch_orth','frac_ue_sizearch_orth','frac_ue_perS_sizearch_orth','macro_ue_sizearch_orth','macro_ue_only_sizearch_orth','pico_ue_sizearch_orth','pico_ue_only_sizearch_orth','Cluster_ue_sizearch_orth','Cluster_ue_only_sizearch_orth','BS_time_sizearch_orth','rate_cvx_sizearch_orth','utility_cvx_sizearch_orth','x_max','act_frac_max','rate_max','utility_max','y_cvx_ic','lambda_cvx_ic','cvx_time_ic','cvx_status_ic','frac_ue_ic','frac_ue_perS_ic','macro_ue_ic','macro_ue_only_ic','pico_ue_ic','pico_ue_only_ic','Cluster_ue_ic','Cluster_ue_only_ic','BS_time_ic','rate_cvx_ic','utility_cvx_ic','x_cvx_BS','rate_cvx_BS','frac_ue_BS','pico_ue_BS','utility_cvx_BS');

end
