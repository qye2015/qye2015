% Function to get the SINR and spectral efficiency 
% 
% ===================================================
% Inputs:
%   BS:     (nB*2 matrix)           BS locations
%   U:      (nU*2 matrix)           user locations
%   L0:     (number)                reference distance to get path loss
%   alpha:  (number)                path loss component
%   M:      (number)                number of BS antennas
%   Uj:     (number)                number of users to be served simutaneously
%   P:      (nB*1 matrix)           transmit power of BSs
% ===================================================
% Outputs:
%   c_zf:   ((2^nB-1)*nU matrix)    spectral efficiency of users using ZF
%   S_zf:   ((2^nB-1)*nU matrix)    received signal using ZF
%   IN_zf:  ((2^nB-1)*nU matrix)    interference + noise using ZF
%   SINR_zf:((2^nB-1)*nU matrix)    SINR using ZF
%   c_mrt:   ((2^nB-1)*nU matrix)   spectral efficiency of users using MRT
%   S_mrt:   ((2^nB-1)*nU matrix)   received signal using MRT
%   IN_mrt:  ((2^nB-1)*nU matrix)   interference + noise using MRT
%   SINR_mrt:((2^nB-1)*nU matrix)   SINR using mrt
%   actBSv: ((2^nB-1)*nB matrix)    the active BSs per cluster
%
% last updated: 6/26/14 2:56pm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [c_zf, S_zf, IN_zf, SINR_zf, c_mrt, S_mrt, IN_mrt, SINR_mrt] = SINR_BS_het(distUB, nM, nP, L0, alpha_m, alpha_p, M, Uj, P)

% =================================================
% initialization
% number of BSs
nB                  =   length(distUB(:, 1));
% number of users
nU                  =   length(distUB(1, :));
% =================================================
% % distance between BSs and users
% distUB              =   pdist2(BS, U, 'euclidean');
% path loss exponent
alpha               =   [alpha_m * ones(nM, 1); alpha_p * ones(nP, 1)];
% path loss matrix: nB*nU
G                   =   repmat(P, 1, nU) ./ (0 + (distUB / L0) .^ repmat(alpha,1,nU));


% =================================================
% ZF
% diversity order factor using ZF
dof                 =   (M - Uj + 1) ./ Uj;

S_zf                =   repmat(dof,1,nU) .* G;
IN_zftp             =  	sum(G, 1) + 1; 
IN_zf               =   repmat(IN_zftp, nB, 1) - G;
SINR_zf             =   S_zf ./ IN_zf;
c_zf                =   log(1 + SINR_zf) / log(2);


% =================================================
% mrt
S_mrt               =   repmat(M ./ Uj, 1, nU) .* G;
IN_mrttp            =   IN_zftp;
IN_mrt              =   (repmat(IN_mrttp, nB, 1) - S_mrt ./ repmat(M, 1, nU));
SINR_mrt            =   S_mrt ./ IN_mrt;
c_mrt               =   log(1 + SINR_mrt) / log(2);


end