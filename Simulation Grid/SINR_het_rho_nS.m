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
% last updated: 12/4/14 3:33pm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [c_zf, S_zf, IN_zf, SINR_zf, c_mrt, S_mrt, IN_mrt, SINR_mrt, Uj_rho, actBSv, archv, n] = SINR_het_rho_nS(NCmax, NCmin, distUB, nM, nP, L0, alpha_m, alpha_p, M, Uj, P, rho)

% =================================================
% initialization
% number of BSs
nB                      =   length(distUB(:, 1));
% number of users
nU                      =   length(distUB(1, :));
% BSs indicator per clusters
yb                      =   zeros(1, nB);
% =================================================
% path loss exponent vector
alpha                   =   [alpha_m * ones(nM, 1); alpha_p * ones(nP, 1)];
% path loss matrix: nB*nU
G                       =   repmat(P, 1, nU) ./ (1 + (distUB / L0) .^ repmat(alpha,1,nU));
IN_zftp                 =  	sum(G, 1) + 1; 
% build a matrix to get the cross term 
cross                   =   zeros(nB, nB, nU);
Uj_rho                  =   ones(nB,1);
for i = 1 : nU
    cross(:,:,i)        =   G(:,i) * G(:,i)';
end

% =================================================
% start clusters    
n                       =   0;
for ntp = 1 : 2^nB - 1
    index               =   str2double(dec2bin(ntp));
    % find the BS indicator
    for j = 1 : nB
        yb(j)           =   mod(round(index/10^(j-1)),2);
    end
    if sum(yb) <= NCmax && sum(yb) >= NCmin
        % number of clusters
        n               =   n + 1;
        % active BS index
        actBS           =   find(yb == 1);
        actBSv(n, actBS)=   1;
        
        % get the arch matrix: which clusters are in which arch
        archv(sum(yb),n)=   1;
        
        % get the signal term
        actmatrix       = zeros(nB, nB);
        actmatrix(actBS, actBS) =   1;

        % diversity order factor using ZF
        Uj_rho(:, n)    =   1;
        Uj_rho(actBS,n) =   round(max(Uj(actBS), Uj(actBS)*rho(length(actBS))));
        dof           	=   (M - Uj_rho(:,n) + 1) ./ Uj_rho(:,n);
        % parameter for mrt
        factor_mrt     	=   M ./ Uj_rho(:,n);
        % =================================================
        % get SINR
        % =================================================
        for i = 1 : nU

            S_zf(n, i) 	=   sum(sum(sqrt(cross(:,:,i) .* (dof * dof')) .* actmatrix));
            S_mrt(n, i)	=    sum(sum(sqrt(cross(:,:,i) .* (factor_mrt * factor_mrt')) .*actmatrix));
            %S_zf(n,i)	=   sum(sum(sqrt(G(:,i) * G(:,i)') * P * dof .* actmatrix));
            %S_mrt(n,i)	= 	sum(sum(sqrt(G(:,i) * G(:,i)') * P * M / Uj .*actmatrix));
        end
        IN_zf(n, :)  	=  	sum(G ,1) - sum(G(actBS, :),1) + 1;
        IN_mrt(n, :)  	=   sum(G, 1) + 1 - sum(1./repmat(Uj_rho(actBS,n),1,nU).*G(actBS,:), 1);

        SINR_zf(n,:)   	=   S_zf(n, :) ./ IN_zf(n, :);
        SINR_mrt(n,:)   =   S_mrt(n, :) ./ IN_mrt(n, :);
        % rate
        c_zf(n,:)       =   log(1 + SINR_zf(n,:)) / log(2);
        c_mrt(n,:)      =   log(1 + SINR_mrt(n, :)) / log(2);
    end

end

end