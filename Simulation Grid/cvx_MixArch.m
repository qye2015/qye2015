% Function to get the cvx for time sharing of archs, 
% where each arch has same-size clus
% 
% ===================================================
% Inputs:
%   NCmax:  (number)        	number of clusters' maximial size
%   NCmin:  (number)            number of clusters' minimal size
%   distUB: (nB*nU matrix)      distance between BSs and users
%   nM:     (number)            number of macro BSs
%   nP:     (number)         	number of pico BSs
%   L0:     (number)        	reference distance to get path loss
%   alpha_m:(number)            path loss component between macro BS & UE
%   alpha_p:(number)            path loss component between pico BS & UE
%   M:      (number)          	number of BS antennas
%   Uj:     (number)           	number of users to be served simutaneously
%   P:      (nB*1 matrix)    	transmit power of BSs
%   rho:    (nC*1 matrix)   	the factor for BSs' capacility of
%   N_nearBS:(number)           the number of nearby BSs that may become a cluster       
%   
% ===================================================
% Outputs:
%   c_zf:   (nC*nU matrix)      spectral efficiency of users using ZF
%   S_zf:   (nC*nU matrix)      received signal using ZF
%   IN_zf:  (nC*nU matrix)      interference + noise using ZF
%   SINR_zf:(nC*nU matrix)      SINR using ZF
%   c_mrt:  (nC*nU matrix)      spectral efficiency of users using MRT
%   S_mrt:  (nC*nU matrix)      received signal using MRT
%   IN_mrt: (nC*nU matrix)      interference + noise using MRT
%   SINR_mrt:(nC*nU matrix)     SINR using mrt
%   Uj_rho: (nB*nC matrix)    	number of users that can be served
%   actBSv: (nC*nB matrix)      active BS indicators
%   archv:  (nA*nC matrix)      cluster indicators in arch
%   n:      (number)            number of clusters
%
% last updated: 1/31/15 11:45am
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [y_cvx_predef, lambda_cvx_predef, cvx_time_predef, cvx_status_predef, frac_ue_predef, frac_ue_perS_predef, macro_ue_predef,  macro_ue_only_predef, pico_ue_predef, pico_ue_only_predef, Cluster_ue_predef, Cluster_ue_only_predef, BS_time_predef, rate_cvx_predef, utility_cvx_predef] = cvx_MixArch(nu, nM, nP, nB, narch, nCluster, SU_js, act_u_cluster, actBSv, archv, c_zf_rho, c_zf_BS)

% =================================================
% some matrixs to be used
% number of clusters per user
nCluster_perU           =   size(act_u_cluster, 2);
% =================================================
% the matrix to get j\in C, |C|=s
actBS_arch              =   zeros(nB*narch, nCluster);
for i = 1 : narch
    actBS_arch(nB*(i-1)+1:nB*i, :) = actBSv'.*repmat(archv(i,:), nB, 1);
end
% =================================================
% the matrix indicating j\in C |C|=s for each user
% col:(j1s1, j2s1, ..., j1s2, j2s2, ...)
% row: c1u1, c2u1, ..., c1u2, c2u2, ...
actBS_arch_u            =   zeros(nB * narch, nCluster_perU * nu);
for i = 1 : nCluster_perU
    actBS_arch_u(:, i:nCluster_perU:nu*nCluster_perU) = actBS_arch(:, act_u_cluster(:, i));
end
% =================================================
% the matrix to get |C|=s
actcluster_arch_u     	=   zeros(narch*nu, nCluster_perU);
for i = 1 : narch
    actcluster_arch_u(nu*(i-1)+1: nu*i, :) = reshape(archv(i, act_u_cluster), nu, nCluster_perU);
end
nCluster_perU_predef    =   nCluster_perU + nB * (narch - 1);
        
% =================================================
% Each predefined arch includes clusters with same size and individual BSs
% rate matrix in the predef arch
c_zf_predef             =   [c_zf_rho; repmat(c_zf_BS, narch-1, 1)];
actBS_arch_u_predef     =   [zeros(nB, nB*(narch-1)*nu); ones(nB, nB*nu), zeros(nB, nB*(narch-2)*nu); zeros(nB,nB*nu), ones(nB,nB*nu), zeros(nB,nB*nu); zeros(nB, nB*(narch-2)*nu), ones(nB, nB*nu)];
actcluster_arch_u_predef=   [zeros(1,nB*(narch-1));ones(1,nB), zeros(1,nB*(narch-2));zeros(1,nB),ones(1,nB),zeros(1,nB);zeros(1,nB*(narch-2)),ones(1,nB)];
% =================================================
% start cvx 
start_time = tic;
cvx_expert true
cvx_clear
cvx_begin
    variable y(nCluster_perU_predef, nu);
    variable lambda(narch, 1);
    maximize( sum(log(sum(c_zf_predef.*y, 1))) );
    %maximize( geo_mean(sum(c_zf_rho.*y, 1)) );
subject to
    actBS_arch_u*reshape(y(1:nCluster_perU,:), nCluster_perU*nu, 1)./SU_js +  actBS_arch_u_predef*reshape(y(1+nCluster_perU:nCluster_perU_predef,:), nB*(narch-1)*nu, 1)./SU_js <= kron(lambda,ones(nB,1));
    sum(actcluster_arch_u.*repmat((y(1:nCluster_perU,:))',narch,1),2) + reshape((actcluster_arch_u_predef*y(nCluster_perU+1:nCluster_perU_predef,:))',nu*narch,1) <= kron(lambda, ones(nu,1));
    %sum(y, 1) <= 1;
    sum(lambda) <= 1;
    y >= 0;
    lambda >= 0;
cvx_end;
cvx_time_predef         =   toc(start_time);
cvx_status_predef       =   cvx_status;

% =================================================
% data post-process
y2                      =   y;
y2(y<10^(-5))           =   0;
y2(y>1-10^(-5))         =   1;
y_cvx_predef           	=   y2;
lambda2              	=   lambda;
lambda2(lambda<10^(-5))	=   0;
lambda2(lambda>1-10^(-5))=   1;
lambda_cvx_predef    	=   lambda2;

% =================================================
% find fractional users -- users connecting to multiple clusters in an arch.
ytp                     =   y_cvx_predef;
ytp(y_cvx_predef>0) 	=   1;
% how many clusters (no individual BS) with each size the user connecting to : nA * nu
ytp_perS_predef_c   	=   (reshape(sum(actcluster_arch_u .* repmat(ytp(1:nCluster_perU,:)', narch, 1),2), nu, narch))';
% how many clusters (only individual BS) with each size the user connecting to : nA * nu
ytp_perS_predef_BS      =   actcluster_arch_u_predef*ytp(nCluster_perU+1:nCluster_perU_predef,:);
% -------------------------------------------
% find the total number of fractional users
frac_ue_id_predef    	=   find(max(max(ytp_perS_predef_c, [], 2),max(ytp_perS_predef_BS, [], 2)) > 1);
frac_ue_predef       	=   length(frac_ue_id_predef);
% find the number of fractional users per s (i.e., clusters with size s)
frac_ue_indic_perS_predef=   (ytp_perS_predef_c > 1 | ytp_perS_predef_BS > 1);
frac_ue_perS_predef    	=   sum(frac_ue_indic_perS_predef, 2);

% =================================================
% matrix to indicate macro BS
macro_index             =   find(sum(actBSv(:,1:nM),2)>0 & sum(actBSv,2) ==1);
actcluster_macro_u      =   zeros(nCluster_perU, nu);
for i = 1 : length(macro_index)
    actcluster_macro_u   = 	min(1, actcluster_macro_u + (act_u_cluster == macro_index(i))');
end
macro_index_pf          =   [];
for i = 1 : (narch-1)
    macro_index_pf      =   [macro_index_pf, nCluster_perU + nB*(i-1)+[1:nM]];
end
% find users only connecting to macro BSs
if ~isempty(macro_index)
    ytp_macro_pf_c      =   ytp(1:nCluster_perU,:) .* actcluster_macro_u;
    macro_uetp          =   find((max(ytp_macro_pf_c,[],1)>0 | sum(ytp(macro_index_pf,:),1)>0) & max(ytp_perS_predef_c, [], 1) <= 1 & max(ytp_perS_predef_BS, [], 1) <= 1);
    macro_ue_predef     =   length(macro_uetp);

    % users that is served by macro BSs in any instance
    macro_ue_only_predef=   length(find((max(ytp_macro_pf_c,[],1)>0 | sum(ytp(macro_index_pf,:),1)>0) & (sum(ytp(1:nCluster_perU,:),1)<=1 | sum(ytp(nCluster_perU+1:nCluster_perU_predef,:),1)<=1)));
end

% =================================================
% matrix to indicate pico BS
pico_index              =   find(sum(actBSv(:,1+nM:nB),2)>0 & sum(actBSv,2) ==1);
actcluster_pico_u       =   zeros(nCluster_perU, nu);
        for i = 1 : length(pico_index)
            actcluster_pico_u   = 	min(1, actcluster_pico_u + (act_u_cluster == pico_index(i))');
        end
pico_index_pf           =   [];
for i = 1 : (narch-1)
    pico_index_pf       =   [pico_index_pf, nCluster_perU + nB*(i-1)+[1+nM:nB]];
end
% find users only connecting to macro BSs
if ~isempty(pico_index)
    ytp_pico_pf_c       =   ytp(1:nCluster_perU,:) .* actcluster_pico_u;
    pico_uetp           =   find((max(ytp_pico_pf_c,[],1)>0 | sum(ytp(pico_index_pf,:),1)>0) & max(ytp_perS_predef_c, [], 1) <= 1 & max(ytp_perS_predef_BS, [], 1) <= 1);
    pico_ue_predef      =   length(pico_uetp);

    % users that is served by macro BSs in any instance
    pico_ue_only_predef =   length(find((max(ytp_pico_pf_c,[],1)>0 | sum(ytp(pico_index_pf,:),1)>0) & (sum(ytp(1:nCluster_perU,:),1)<=1 | sum(ytp(nCluster_perU+1:nCluster_perU_predef,:),1)<=1)));
end  

% =================================================
% matrix to indicate cluster
cluster_index          	=   find(sum(actBSv,2) > 1);
actcluster_cluster_u 	=   zeros(nCluster_perU, nu);
        for i = 1 : length(cluster_index)
            actcluster_cluster_u=   min(1, actcluster_cluster_u + (act_u_cluster == cluster_index(i))');
        end
% find users only connecting to macro BSs
if ~isempty(cluster_index)
    ytp_cluster_pf_c  	=   ytp(1:nCluster_perU,:) .* actcluster_cluster_u;
    Cluster_uetp      	=   find(max(ytp_cluster_pf_c,[],1)>0  & max(ytp_perS_predef_c, [], 1) <= 1 & max(ytp_perS_predef_BS, [], 1) <= 1);
    Cluster_ue_predef 	=   length(Cluster_uetp);

    % users that is served by macro BSs in any instance
    Cluster_ue_only_predef=	length(find(max(ytp_cluster_pf_c,[],1)>0 & sum(ytp(1:nCluster_perU,:),1)<=1 ));
end

% =================================================
% find how long each BS working in clusters
BS_time_predef          =   reshape(actBS_arch_u*reshape(y_cvx_predef(1:nCluster_perU, :), nCluster_perU*nu, 1)./SU_js, nB, narch);
for i = 1 : nB
    indextp             =   nCluster_perU + [i:nB:nB*(narch-1)];
    BS_time_predef(i,1) =   BS_time_predef(i, 1) + sum(sum(y_cvx_predef(indextp,:)))/SU_js(i);
end

% =================================================
% utility function
rate_cvx_predef         =   sum(c_zf_predef.*y_cvx_predef, 1);
utility_cvx_predef  	=   sum(log(rate_cvx_predef));

end