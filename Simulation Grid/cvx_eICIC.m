% Function to get the cvx results for eICIC
% eICIC: two architectures - cellular mode and cellular mode with off macro
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
function [y_cvx_ic, lambda_cvx_ic, cvx_time_ic, cvx_status_ic, frac_ue_ic, frac_ue_perS_ic, macro_ue_ic,  macro_ue_only_ic, pico_ue_ic, pico_ue_only_ic, Cluster_ue_ic, Cluster_ue_only_ic, BS_time_ic, rate_cvx_ic, utility_cvx_ic] = cvx_eICIC(nu, nM, nP, nB, nCluster, narch, SU_js, SU, act_u_cluster, actBSv, archv, c_zf_rho, c_zf_off)
% ===================================================
% some matrixs to be used
% number of clusters per user
nCluster_perU                   =   size(act_u_cluster, 2);
nCluster_perU_ic                =   nCluster_perU + nP;
% number of clusters 
nCluster_ic                     =   nCluster + nP;
% number of architectures: narch + off arch
narch_ic                        =   narch + 1;
% ===================================================
% the vector of S_{js}: (nB*narch)*1
SU_ic                           =   [SU_js; SU];
% ===================================================
% spectral efficiency
c_zf_ic                         =   [c_zf_rho; c_zf_off(nM+1:nB, :)];
% ===================================================
% the matrix to get j\in C, C\in A
actBS_arch                      =   zeros(nB*narch, nCluster);
for i = 1 : narch
    actBS_arch(nB*(i-1)+1:nB*i, :) = actBSv'.*repmat(archv(i,:), nB, 1);
end
actBS_arch_ic                   =   zeros(nB*narch_ic, nCluster_ic);
actBS_arch_ic(1:nB*narch, 1:nCluster)=  actBS_arch;
actBS_arch_ic(nB*narch+nM+1 : nB*narch_ic, nCluster+1:nCluster_ic) =   eye(nP);
% ===================================================
% matrix indicating which user is served by which cluster
act_u_cluster_ic                =   zeros(nu, nCluster_perU_ic);
act_u_cluster_ic(:,1:nCluster_perU) =   act_u_cluster;
act_u_cluster_ic(:,nCluster_perU+1:nCluster_perU_ic) =  ones(nu, nP);
% =================================================
% the matrix indicating j\in C |C|=s for each user
% col:(j1s1, j2s1, ..., j1s2, j2s2, ...)
% row: c1u1, c2u1, ..., c1u2, c2u2, ...
% actBS_arch_u            =   zeros(nB * narch, nCluster_perU * nu);
% for i = 1 : nCluster_perU
%     actBS_arch_u(:, i:nCluster_perU:nu*nCluster_perU) = actBS_arch(:, act_u_cluster(:, i));
% end
actBS_arch_u_ic         =   zeros(nB * narch_ic, nCluster_perU_ic * nu);
for i = 1 : nCluster_perU_ic
    actBS_arch_u_ic(:, i:nCluster_perU_ic:nu*nCluster_perU_ic) = actBS_arch_ic(:, act_u_cluster_ic(:, i));
end
% =================================================
% the matrix to get C\in A
actcluster_arch_u     	=   zeros(narch*nu, nCluster_perU);
for i = 1 : narch
    actcluster_arch_u(nu*(i-1)+1: nu*i, :) = reshape(archv(i, act_u_cluster), nu, nCluster_perU);
end  
actcluster_arch_u_ic    =   zeros(narch_ic*nu, nCluster_perU_ic);
actcluster_arch_u_ic(1:narch*nu, 1:nCluster_perU) = actcluster_arch_u;
actcluster_arch_u_ic(narch*nu+1:narch_ic*nu, nCluster_perU+1:nCluster_perU_ic) =  ones(nu, nP);

             
% =================================================
% Each predefined arch includes clusters with same size 
% start cvx 
start_time = tic;
cvx_expert true
cvx_clear
cvx_begin
    variable y(nCluster_perU_ic, nu);
    variable lambda(narch_ic, 1);
    maximize( sum(log(sum(c_zf_ic.*y, 1))) );
subject to
    actBS_arch_u_ic*reshape(y, nCluster_perU_ic*nu, 1)./SU_ic <= kron(lambda,ones(nB,1));
    sum(actcluster_arch_u_ic.*repmat(y',narch_ic,1),2) <= kron(lambda, ones(nu,1));
    %sum(y, 1) <= 1;
   	sum(lambda) <= 1;
  	y >= 0;
  	lambda >= 0;
cvx_end;
cvx_time_ic             =   toc(start_time);
cvx_status_ic           =   cvx_status;
% =================================================
% data post-process
y2                      =   y;
y2(y<10^(-5))           =   0;
y2(y>1-10^(-5))         =   1;
y_cvx_ic                =   y2;
lambda2              	=   lambda;
lambda2(lambda<10^(-5))	=   0;
lambda2(lambda>1-10^(-5))=   1;
lambda_cvx_ic           =   lambda2;
        
% =================================================
% find fractional users -- users connecting to multiple clusters in an arch.
ytp                     =   y_cvx_ic;
ytp(y_cvx_ic>0)         =   1;
% how many clusters (no individual BS) with each size the user connecting to : nA * nu
ytp_perS_ic             =   (reshape(sum(actcluster_arch_u_ic .* repmat(ytp', narch_ic, 1),2), nu, narch_ic))';
% -------------------------------------------
% find the total number of fractional users
frac_ue_id_ic           =   find(max(ytp_perS_ic, [], 2) > 1);
frac_ue_ic              =   length(frac_ue_id_ic);
% find the number of fractional users per s (i.e., clusters with
% size s)
frac_ue_indic_perS_ic   =   ytp_perS_ic > 1;
frac_ue_perS_ic         =   sum(frac_ue_indic_perS_ic, 2);
        
% =================================================
% find users only connecting to macro BSs
macro_index             =   find(sum(actBSv(:,1:nM),2)>0 & sum(actBSv,2) ==1);
actcluster_macro_u      =   zeros(nCluster_perU_ic, nu);
for i = 1 : length(macro_index)
    actcluster_macro_u   = 	min(1, actcluster_macro_u + (act_u_cluster_ic == macro_index(i))');
end
if ~isempty(macro_index)
    ytp_macro_ic        =   ytp .* actcluster_macro_u;
    macro_uetp          =   find(max(ytp_macro_ic,[],1)>0  & max(ytp_perS_ic, [], 1) <= 1);
    macro_ue_ic         =   length(macro_uetp);
    % users that is served by macro BSs in any instance
    macro_ue_only_ic    =   length(find(max(ytp_macro_ic,[],1)>0 & sum(ytp,1)<=1 ));
end
        
% =================================================
% find users only connecting to macro BSs
pico_index              =   find(sum(actBSv(:,1+nM:nB),2)>0 & sum(actBSv,2) ==1);
actcluster_pico_u       =   zeros(nCluster_perU_ic, nu);
for i = 1 : length(pico_index)
    actcluster_pico_u   = 	min(1, actcluster_pico_u + (act_u_cluster_ic == pico_index(i))');
end
if ~isempty(pico_index)
    ytp_pico_ic         =   ytp .* actcluster_pico_u;
    pico_uetp           =   find(max(ytp_pico_ic,[],1)>0 & max(ytp_perS_ic, [], 1) <= 1);
    pico_ue_ic          =   length(pico_uetp);
    % users that is served by macro BSs in any instance
    pico_ue_only_ic     =   length(find(max(ytp_pico_ic,[],1)>0 & sum(ytp,1)<=1 ));
end  
        
% =================================================
% find users only connecting to macro BSs
cluster_index          	=   find(sum(actBSv,2) > 1);
actcluster_cluster_u 	=   zeros(nCluster_perU_ic, nu);
for i = 1 : length(cluster_index)
    actcluster_cluster_u=   min(1, actcluster_cluster_u + (act_u_cluster_ic == cluster_index(i))');
end
if ~isempty(cluster_index)
    ytp_cluster_ic      =   ytp .* actcluster_cluster_u;
    Cluster_uetp      	=   find(max(ytp_cluster_ic,[],1)>0  & max(ytp_perS_ic, [], 1) <= 1);
    Cluster_ue_ic       =   length(Cluster_uetp);
    % users that is served by macro BSs in any instance
    Cluster_ue_only_ic  =	length(find(max(ytp_cluster_ic,[],1)>0 & sum(ytp,1)<=1 ));
end
        
% =================================================
% find how long each BS working in clusters
BS_time_ic              =   reshape(actBS_arch_u_ic*reshape(y_cvx_ic, nCluster_perU_ic*nu, 1)./SU_ic, nB, narch_ic);

% =================================================
% utility function
rate_cvx_ic             =   sum(c_zf_ic.*y_cvx_ic, 1);
utility_cvx_ic          =   sum(log(rate_cvx_ic));

   
            
end