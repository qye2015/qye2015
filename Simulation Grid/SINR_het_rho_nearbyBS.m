% Function to get the SINR and spectral efficiency
% Cluster with size nmin-nmax
% Each architecture includes clusters with the same size
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
% last updated: 1/29/15 11:45am
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [c_zf, S_zf, IN_zf, SINR_zf, c_mrt, S_mrt, IN_mrt, SINR_mrt, Uj_rho, actBSv, archv, act_u_cluster, ncluster] = SINR_het_rho_nearbyBS(NCmax, NCmin, distUB, nM, nP, L0, alpha_m, alpha_p, M, Uj, P, rho, N_nearBS, c_zf_BS)

% =================================================
% initialization
% number of BSs
nB                      =   length(distUB(:, 1));
% number of users
nu                      =   length(distUB(1, :));
% BSs indicator per clusters
yb                      =   zeros(1, nB);
% =================================================
% path loss exponent vector
alpha                   =   [alpha_m * ones(nM, 1); alpha_p * ones(nP, 1)];
% path loss matrix: nB*nU
G                       =   repmat(P, 1, nu) ./ (1 + (distUB / L0) .^ repmat(alpha,1,nu));
IN_zftp                 =  	sum(G, 1) + 1; 
% build a matrix to get the cross term 
cross                   =   zeros(nB, nB, nu);
Uj_rho                  =   ones(nB,1);
for i = 1 : nu
    cross(:,:,i)        =   G(:,i) * G(:,i)';
end

% =================================================
% the number of clusters to be considered per UE
N_cluster_perU          =   0;
for i = NCmin : NCmax
    N_cluster_perU    	=   N_cluster_perU + factorial(N_nearBS)/factorial(i)/factorial(N_nearBS-i);
end

act_u_cluster           =   zeros(nu, N_cluster_perU);
actBSv                  =   zeros(1, nB);
archv                   =   zeros(NCmax, 1);
ncluster                =   0;
% =================================================
% start clusters for each user
for i                   =   1 : nu
    nctp                =   0;
    % find the nearest N_nearBS BSs
    [~, sort_index]     =   sort(c_zf_BS(:, i), 'descend');
    BS_index            =   sort_index(1 : N_nearBS);
    
    % ----------------------------------------------
    % find the clusters per user
    for ntp = 1 : 2^N_nearBS - 1
        index                   =   str2double(dec2bin(ntp));
        % find the BS indicator
        yb                      =   zeros(1, nB);
        for j = 1 : N_nearBS
            yb(BS_index(j))  	=   mod(round(index/10^(j-1)),2);
        end
        % ----------------------------------------------
        % check if the size of clusters is smaller than NCmax
        if sum(yb) <= NCmax && sum(yb) >= NCmin
            % the BS index
            actBS_index         =   find(yb == 1);
            % check whether this cluster has been recorded or not
            if min(sum(abs(actBSv - repmat(yb, length(actBSv(:,1)),1)),2)) > 0 % not recorded
                ncluster      	=    ncluster + 1;
                cluster_index 	=    ncluster;
                actBSv(cluster_index, actBS_index)  =	1;
            end
            %cluster_index    	=   sum(2.^(find(yb==1)-1));
            cluster_index       =   find(sum(abs(actBSv - repmat(yb, length(actBSv(:,1)),1)),2)==0);
            nctp                =   nctp + 1;
            act_u_cluster(i, nctp)  =   cluster_index;
            
            actBSv(cluster_index, actBS_index)      =   1;
            archv(sum(yb), cluster_index) = 1;
            
            % ----------------------------------------------
            % get the signal term
            actmatrix           = zeros(nB, nB);
            actmatrix(actBS_index, actBS_index) =   1;

            % diversity order factor using ZF
            Uj_rho(:, cluster_index)    =   1;
            Uj_rho(actBS_index,cluster_index) =   round(max(Uj(actBS_index), Uj(actBS_index)*rho(length(actBS_index))));
            dof                 =   (M - Uj_rho(:, cluster_index) + 1) ./ Uj_rho(:,cluster_index);
            % parameter for mrt
            factor_mrt          =   M ./ Uj_rho(:, cluster_index);
            

            S_zf(nctp, i)       =   sum(sum(sqrt(cross(:,:,i) .* (dof * dof')) .* actmatrix));
            S_mrt(nctp, i)      =    sum(sum(sqrt(cross(:,:,i) .* (factor_mrt * factor_mrt')) .*actmatrix));
            %S_zf(n,i)	=   sum(sum(sqrt(G(:,i) * G(:,i)') * P * dof .* actmatrix));
            %S_mrt(n,i)	= 	sum(sum(sqrt(G(:,i) * G(:,i)') * P * M / Uj .*actmatrix));
        
            IN_zf(nctp, i)  	=  	sum(G(:, i)) - sum(G(actBS_index, i)) + 1;
            IN_mrt(nctp, i) 	=   sum(G(:, i)) + 1 - sum(1./(Uj_rho(actBS_index,cluster_index)).*G(actBS_index,i), 1);

            SINR_zf(nctp, i)  	=   S_zf(nctp, i) ./ IN_zf(nctp, i);
            SINR_mrt(nctp, i) 	=   S_mrt(nctp, i) ./ IN_mrt(nctp, i);
            % rate
            c_zf(nctp, i)   	=   log(1 + SINR_zf(nctp, i)) / log(2);
            c_mrt(nctp, i)  	=   log(1 + SINR_mrt(nctp, i)) / log(2);
        end
    end
    
    % output
    if mod(i,100) == 0
        fprintf('finish calculation for nu=%d, ncluster=%d\n', i,ncluster);
    end
    
end


% % =================================================
% % start clusters    
% n                       =   0;
% for ntp = 1 : 2^nB - 1
%     index               =   str2double(dec2bin(ntp));
%     % find the BS indicator
%     for j = 1 : nB
%         yb(j)           =   mod(round(index/10^(j-1)),2);
%     end
%     if sum(yb) <= NCmax && sum(yb) >= NCmin
%         % number of clusters
%         n               =   n + 1;
%         % active BS index
%         actBS           =   find(yb == 1);
%         actBSv(n, actBS)=   1;
%         
%         % get the arch matrix: which clusters are in which arch
%         archv(sum(yb),n)=   1;
%         
%         % get the signal term
%         actmatrix       = zeros(nB, nB);
%         actmatrix(actBS, actBS) =   1;
% 
%         % diversity order factor using ZF
%         Uj_rho(:, n)    =   1;
%         Uj_rho(actBS,n) =   round(max(Uj(actBS), Uj(actBS)*rho(length(actBS))));
%         dof           	=   (M - Uj_rho(:,n) + 1) ./ Uj_rho(:,n);
%         % parameter for mrt
%         factor_mrt     	=   M ./ Uj_rho(:,n);
%         % =================================================
%         % get SINR
%         % =================================================
%         for i = 1 : nU
% 
%             S_zf(n, i) 	=   sum(sum(sqrt(cross(:,:,i) .* (dof * dof')) .* actmatrix));
%             S_mrt(n, i)	=    sum(sum(sqrt(cross(:,:,i) .* (factor_mrt * factor_mrt')) .*actmatrix));
%             %S_zf(n,i)	=   sum(sum(sqrt(G(:,i) * G(:,i)') * P * dof .* actmatrix));
%             %S_mrt(n,i)	= 	sum(sum(sqrt(G(:,i) * G(:,i)') * P * M / Uj .*actmatrix));
%         end
%         IN_zf(n, :)  	=  	sum(G ,1) - sum(G(actBS, :),1) + 1;
%         IN_mrt(n, :)  	=   sum(G, 1) + 1 - sum(1./repmat(Uj_rho(actBS,n),1,nU).*G(actBS,:), 1);
% 
%         SINR_zf(n,:)   	=   S_zf(n, :) ./ IN_zf(n, :);
%         SINR_mrt(n,:)   =   S_mrt(n, :) ./ IN_mrt(n, :);
%         % rate
%         c_zf(n,:)       =   log(1 + SINR_zf(n,:)) / log(2);
%         c_mrt(n,:)      =   log(1 + SINR_mrt(n, :)) / log(2);
%     end
% 
% end

end