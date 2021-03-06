% Function to get the cvx results for BS arch,
% where clusters with different sizes can work together
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
function [y_cvx_mixA, mu_cvx_mixA, cvx_time_mixA, cvx_status_mixA, frac_ue_mixA, frac_ue_perS_mixA, macro_ue_mixA,  macro_ue_only_mixA, pico_ue_mixA, pico_ue_only_mixA, cluster_ue_mixA, cluster_ue_only_mixA, BS_time_mixA, rate_cvx_mixA, utility_cvx_mixA] = cvx_mixBSarch_orth(nu, nB, narch, nCluster, SU_js, act_u_cluster, actBSv, archv, c_zf_rho)

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
% =================================================
% matrix to indicator 
% row: j1s1(u1, u2, u3, ...), j2s1(u1,u2,u3,...), ..., j1s2(u1,u2,u3...), j2s2(u1,u2,u3)
% col: c1, c2, c3, ...
act_arch_u_cluster      =   reshape((reshape(actBS_arch_u,[],nu))',[], nCluster_perU);
% =================================================
% start cvx 
start_time = tic;
cvx_expert true
cvx_clear
cvx_begin
    variable y(nCluster_perU, nu);
    variable mu(nB, narch);
    maximize( sum(log(sum(c_zf_rho.*y, 1))) );
    %maximize( geo_mean(sum(c_zf_rho.*y, 1)) );
subject to
    actBS_arch_u*reshape(y, nCluster_perU*nu, 1)./SU_js <= reshape(mu,nB*narch,1);
    sum(act_arch_u_cluster.*repmat(y',nB*narch,1),2) <= kron(reshape(mu,nB*narch,1), ones(nu,1));
    sum(y, 1) <= 1;
    sum(mu, 2) <= 1;
    y >= 0;
    mu >= 0;
cvx_end;
cvx_time_mixA           =   toc(start_time);
cvx_status_mixA         =   cvx_status;

% =================================================
% data post-process
y2                      =   y;
y2(y<10^(-5))           =   0;
y2(y>1-10^(-5))         =   1;
y_cvx_mixA              =   y2;
mu2                     =   mu;
mu2(mu<10^(-5))         =   0;
mu2(mu>1-10^(-5))       =   1;
mu_cvx_mixA             =   mu2;

% =================================================
% find fractional users -- users connecting to multiple clusters in an arch.
ytp                     =   y_cvx_mixA;
ytp(y_cvx_mixA>0)       =   1;
% how many clusters with each size the user connecting to: nA*nu
ytp_perS                =   (reshape(sum(actcluster_arch_u .* repmat(ytp', narch, 1),2), nu, narch))';
% -------------------------------------------
% find the total number of fractional users
frac_ue_id              =   find(max(ytp_perS, [], 2) > 1);
frac_ue_mixA            =   length(frac_ue_id);
% find the number of fractional users per s (i.e., clusters with size s)
frac_ue_indic_perS      =   ytp_perS > 1;
frac_ue_perS_mixA       =   sum(frac_ue_indic_perS, 2);
% =================================================
% matrix to indicate macro BS
macro_index             =   find(sum(actBSv(:,1:nM),2)>0 & sum(actBSv,2) ==1);
actcluster_macro_u      =   zeros(nCluster_perU, nu);
for i = 1 : length(macro_index)
    actcluster_macro_u   = 	min(1, actcluster_macro_u + (act_u_cluster == macro_index(i))');
end
% find users only connecting to macro BSs
if ~isempty(macro_index)
    ytp_macro           =   ytp .* actcluster_macro_u;
    macro_uetp          =   find(max(ytp_macro,[],1)>0 & max(ytp_perS, [], 1) <= 1);
    macro_ue_mixA       =   length(macro_uetp);

    % users that is served by macro BSs in any instance
    macro_ue_only_mixA	=   length(find(max(ytp_macro,[],1)>0 & sum(ytp,1)<=1));
end
% =================================================
% matrix to indicate pico BS
pico_index              =   find(sum(actBSv(:,1+nM:nB),2)>0 & sum(actBSv,2) ==1);
actcluster_pico_u       =   zeros(nCluster_perU, nu);
for i = 1 : length(pico_index)
    actcluster_pico_u   = 	min(1, actcluster_pico_u + (act_u_cluster == pico_index(i))');
end
% find users only connecting to small BSs
if ~isempty(pico_index)
    ytp_pico            =   ytp .* actcluster_pico_u;
    pico_uetp           =   find(max(ytp_pico,[],1)>0 & max(ytp_perS, [], 1) <= 1);
    pico_ue_mixA        =   length(pico_uetp);

    % users that is served by macro BSs in any instance
    pico_ue_only_mixA	=   length(find(max(ytp_pico,[],1)>0 & sum(ytp,1)<=1));
end
% =================================================
% matrix to indicate clusters
cluster_index          	=   find(sum(actBSv,2) > 1);
actcluster_cluster_u 	=   zeros(nCluster_perU, nu);
for i = 1 : length(cluster_index)
    actcluster_cluster_u=   min(1, actcluster_cluster_u + (act_u_cluster == cluster_index(i))');
end
% find users only connecting to clusters with multiple BSs
if ~isempty(cluster_index)
    ytp_cluster       	=   ytp .* actcluster_cluster_u;
    cluster_uetp     	=   find(max(ytp_cluster,[],1)>0 & max(ytp_perS, [], 1) <= 1);
    cluster_ue_mixA    	=   length(cluster_uetp);

    % users that is served by macro BSs in any instance
    cluster_ue_only_mixA=   length(find(max(ytp_cluster,[],1)>0 & sum(ytp,1)<=1));
end

% =================================================
% find how long each BS working in clusters
BS_time_mixA            =   mu_cvx_mixA;

% =================================================
% utility function
rate_cvx_mixA           =   sum(c_zf_rho.*y_cvx_mixA, 1);
utility_cvx_mixA        =   sum(log(rate_cvx_mixA));


end