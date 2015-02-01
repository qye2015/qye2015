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
function [y_cvx_sizearch, lambda_cvx_sizearch, cvx_time_sizearch, cvx_status_sizearch, frac_ue_sizearch, frac_ue_perS_sizearch, macro_ue_sizearch,  macro_ue_only_sizearch, pico_ue_sizearch, pico_ue_only_sizearch, Cluster_ue_sizearch, Cluster_ue_only_sizearch, BS_time_sizearch, rate_cvx_sizearch, utility_cvx_sizearch] = cvx_SamesizeArch(nu, nM, nP, nB, narch, nCluster, SU_js, act_u_cluster, actBSv, archv, c_zf_rho)

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
        % Each predefined arch includes clusters with same size 
        % start cvx 
        start_time = tic;
        cvx_expert true
        cvx_clear
        cvx_begin
            variable y(nCluster_perU, nu);
            variable lambda(narch, 1);
            maximize( sum(log(sum(c_zf_rho.*y, 1))) );
            %maximize( geo_mean(sum(c_zf_rho.*y, 1)) );
        subject to
            actBS_arch_u*reshape(y, nCluster_perU*nu, 1)./SU_js <= kron(lambda,ones(nB,1));
            sum(actcluster_arch_u.*repmat(y',narch,1),2) <= kron(lambda, ones(nu,1));
            %sum(y, 1) <= 1;
            sum(lambda) <= 1;
            y >= 0;
            lambda >= 0;
        cvx_end;
        cvx_time_sizearch 	=   toc(start_time);
        cvx_status_sizearch =   cvx_status;
        % =================================================
        % data post-process
        y2                      =   y;
        y2(y<10^(-5))           =   0;
        y2(y>1-10^(-5))         =   1;
        y_cvx_sizearch       	=   y2;
        lambda2              	=   lambda;
        lambda2(lambda<10^(-5))	=   0;
        lambda2(lambda>1-10^(-5))=   1;
        lambda_cvx_sizearch   	=   lambda2;
        
        % =================================================
        % find fractional users -- users connecting to multiple clusters in an arch.
        ytp                     =   y_cvx_sizearch;
        ytp(y_cvx_sizearch>0) 	=   1;
        % how many clusters (no individual BS) with each size the user connecting to : nA * nu
        ytp_perS_sizearch   	=   (reshape(sum(actcluster_arch_u .* repmat(ytp', narch, 1),2), nu, narch))';
        % -------------------------------------------
        % find the total number of fractional users
        frac_ue_id_sizearch    	=   find(max(ytp_perS_sizearch, [], 2) > 1);
        frac_ue_sizearch       	=   length(frac_ue_id_sizearch);
        % find the number of fractional users per s (i.e., clusters with
        % size s)
        frac_ue_indic_perS_sizearch=   ytp_perS_sizearch > 1;
        frac_ue_perS_sizearch  	=   sum(frac_ue_indic_perS_sizearch, 2);
        
        % =================================================
        % find users only connecting to macro BSs
        macro_index             =   find(sum(actBSv(:,1:nM),2)>0 & sum(actBSv,2) ==1);
        actcluster_macro_u      =   zeros(nCluster_perU, nu);
        for i = 1 : length(macro_index)
            actcluster_macro_u   = 	min(1, actcluster_macro_u + (act_u_cluster == macro_index(i))');
        end
        if ~isempty(macro_index)
            ytp_macro_sizearch 	=   ytp .* actcluster_macro_u;
            macro_uetp          =   find(max(ytp_macro_sizearch,[],1)>0  & max(ytp_perS_sizearch, [], 1) <= 1);
            macro_ue_sizearch  	=   length(macro_uetp);
            % users that is served by macro BSs in any instance
            macro_ue_only_sizearch=   length(find(max(ytp_macro_sizearch,[],1)>0 & sum(ytp,1)<=1 ));
        end
        
        % =================================================
        % find users only connecting to macro BSs
        pico_index              =   find(sum(actBSv(:,1+nM:nB),2)>0 & sum(actBSv,2) ==1);
        actcluster_pico_u       =   zeros(nCluster_perU, nu);
        for i = 1 : length(pico_index)
            actcluster_pico_u   = 	min(1, actcluster_pico_u + (act_u_cluster == pico_index(i))');
        end
        if ~isempty(pico_index)
            ytp_pico_sizearch  	=   ytp .* actcluster_pico_u;
            pico_uetp           =   find(max(ytp_pico_sizearch,[],1)>0 & max(ytp_perS_sizearch, [], 1) <= 1);
            pico_ue_sizearch   	=   length(pico_uetp);

            % users that is served by macro BSs in any instance
            pico_ue_only_sizearch =   length(find(max(ytp_pico_sizearch,[],1)>0 & sum(ytp,1)<=1 ));
        end  
        
        % =================================================
        % find users only connecting to macro BSs
        cluster_index          	=   find(sum(actBSv,2) > 1);
        actcluster_cluster_u 	=   zeros(nCluster_perU, nu);
        for i = 1 : length(cluster_index)
            actcluster_cluster_u=   min(1, actcluster_cluster_u + (act_u_cluster == cluster_index(i))');
        end
        if ~isempty(cluster_index)
            ytp_cluster_sizearch=   ytp .* actcluster_cluster_u;
            Cluster_uetp      	=   find(max(ytp_cluster_sizearch,[],1)>0  & max(ytp_perS_sizearch, [], 1) <= 1);
            Cluster_ue_sizearch =   length(Cluster_uetp);
            % users that is served by macro BSs in any instance
            Cluster_ue_only_sizearch=	length(find(max(ytp_cluster_sizearch,[],1)>0 & sum(ytp,1)<=1 ));
        end
        
        % =================================================
        % find how long each BS working in clusters
        BS_time_sizearch     	=   reshape(actBS_arch_u*reshape(y_cvx_sizearch, nCluster_perU*nu, 1)./SU_js, nB, narch);

        % =================================================
        % utility function
        rate_cvx_sizearch       =   sum(c_zf_rho.*y_cvx_sizearch, 1);
        utility_cvx_sizearch	=   sum(log(rate_cvx_sizearch));



end