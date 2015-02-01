%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% get spectral efficiency for all realizations
% Last update 1/27/15 3:00pm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; clf; close all;
% ------------------------------------------
step_rho                            =   1; % 0.75 0.5 0.25];  
% ------------------------------------------
load(strcat('data/',int2str(1),'-nphot-',int2str(3),'-nuhot-',int2str(8),'-grid.mat'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% run data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
loop                =   0;
save('data/loop.mat','loop');
    
while loop < totloop
        
        load('data/loop.mat');
        loop                        =   loop+1;
        if mod(loop, 4) ==0
            clc;
        end
        fprintf('start loop=%d\n',loop);
        
        load(strcat('data/',int2str(loop),'-nphot-',int2str(denseP_hot/denseP_reg),'-nuhot-',int2str(denseu_hot/denseu_reg),'-grid.mat'));
        

        % =================================================
        % MIMO parameters                       
        % user factor
        rho                         =   ones(nB, 1);
        for i = 2 : nB
            rho(i)                  =   step_rho * i;
        end
        
        % =================================================
        % transmit power
        Pt                          =   zeros(nM + nP, 1);
        Pt(1: nM)                   =   P_macro;
        Pt(nM+1 : nM+nP)            =   P_pico;
        % transmit power in eICIC blank RBs
        Pt_ic                       =   zeros(nM + nP, 1);
        %Pt_ic(1 : nM)               =   0;
        Pt_ic(nM+1 : nM+nP)         =   P_pico;
        % =================================================
        % number of antennas
        M                           =   ones(nM + nP, 1);
        M(1 : nM)                   =   M_macro;
        M(nM+1 : nM+nP)             =   M_pico;
        % =================================================
        % number of users that can be served
        SU                          =   ones(nM + nP, 1);
        SU(1 : nM)                  =   S_macro;
        SU(nM+1 : nM+nP)            =   S_pico;
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Get SINR from all possible clusters, with maximial cluster size NC_max
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ===================================================
        % Inputs:
        %   NCmax:  (number)        	number of clusters' maximial size
        %   NCmin:  (number)            number of clusters' minimal size
        %   distBU: (nB*nU matrix)      distance between BSs and users
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
        % =================================================
        
        % get SINR, where each BS in cluster with size n can serve nS users
        [c_zf_rho, ~, ~, SINR_zf_rho, c_mrt_rho, S_mrt, IN_mrt, SINR_mrt, SU_rho, actBSv, archvtp, nCluster] = SINR_het_rho_nS(NC_max, NC_min, DistBU, nM, nP, L0, alpha_macro, alpha_pico, M, SU, Pt, rho);
        fprintf('finish peak rate,  nCluster=%d\n', nCluster);


        % =================================================
        % find the architecture-cluster indicator matrix
        archv                     =   archvtp(NC_min:NC_max,:);
        % number of architectures
        narch                     =   size(archv, 1);


        % =================================================
        % find the index of macro BSs, pico BSs, and clusters with multiple BSs
%         macro_index               =   find(sum(actBSv(:,1:nM),2)==1 & sum(actBSv,2)==1);
%         pico_index                =   find(sum(actBSv(:,nM+1:nB),2)==1 & sum(actBSv,2)==1);
%         BS_index                  =   [macro_index; pico_index];
%         cluster_index             =   find(sum(actBSv,2)>1);


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Get SINR when the macro BSs and pico are using orthogonal bands
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get SINR from picos
        [c_zf_orthtp, ~, ~, SINR_zf_orth, c_mrt_orth, S_mrt_orth, IN_mrt_orth, SINR_mrt_orth, SU_rho_orthtp, actBSv_orthtp, archv_orthtp, nCluster_orthtp] = SINR_het_rho_nS(NC_max, NC_min, DistBU(nM+1:nB,:), 0, nP, L0, alpha_macro, alpha_pico, M(nM+1:nB), SU(nM+1:nB), Pt(nM+1:nB), rho);
        
        % =================================================
        % get SINR from macro BSs, where macro BSs does not become clusters
        [c_zf_macro, S_zf_macro, IN_zf_macro, SINR_zf_macro, c_mrt_macro, S_mrt_macro, IN_mrt_macro, SINR_mrt_macro] = SINR_BS_het(DistBU(1:nM, :), nM, 0, L0, alpha_macro, alpha_pico, M(1:nM), SU(1:nM), Pt(1:nM));
        
        c_zf_orth                   =   [c_zf_macro; c_zf_orthtp];
        
        fprintf('finish peak rate in orthogonal case, nCluster=%d\n', nCluster_orthtp+nM);

        narch_orth                  =   length(archv_orthtp) + 1;
        nCluster_orth               =   nCluster_orthtp + nM;
        SU_rho_orth                 =   zeros(nB, nCluster_orth);
        SU_rho_orth(1:nM,1:nM)      =   SU_rho(1:nM, macro_index);
        SU_rho_orth(nM+1:nB,nM+1:nCluster_orth) =   SU_rho_orthtp;

        actBSv_orth                 =   zeros(nCluster_orth, nB);
        actBSv_orth(1:nM,1:nM)      =   actBSv(macro_index, 1:nM);
        actBSv_orth(nM+1:nCluster_orth,nM+1:nB) = actBSv_orthtp;

        archv_orth                  =   zeros(narch_orth, nCluster_orthtp);
        archv_orth(1,1:nM)          =   archv(1, macro_index);
        archv_orth(2:narch_orth,nM+1:nCluster_orth) = archv_orthtp(NC_min:NC_max,:);

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Get SINR when the macro BSs are off
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % =================================================
        % spectral efficiency when macro BSs are off
        [c_zf_off, S_zf_off, IN_zf_off, SINR_zf_off, c_mrt_off, S_mrt_off, IN_mrt_off, SINR_mrt_off] = SINR_BS_het(DistBU, nM, nP, L0, alpha_macro, alpha_pico, M, SU, Pt_ic);
        fprintf('finish peak rate in eICIC case for static users\n');

        % spectral efficiency in arch2+eICIC
        %c_zf_all_ic         	=   [c_zf_rho;c_zf_off];

    
%             % =================================================
%             % find the architecture-cluster indicator matrix in arch2+eICIC
%             % we let pico BSs in the off RBs become individual clusters
%             archv_all_ic            =   zeros(narch+1, nCluster+nB);
%             archv_all_ic(1:narch,1:nCluster) =	archv;
%             archv_all_ic(narch+1,nCluster+2:nCluster+nB) =  archv(1,pico_index); 
%             % number of architectures
%             narch_all_ic         	=   size(archv_all_ic, 1);
%             % number of clusters in arch2 + eICIC
%             nCluster_all_ic         =   size(archv_all_ic, 2);
%             % indicator of which BS is active in which cluster
%             actBSv_all_ic           =   zeros(nCluster_all_ic, nB);
%             actBSv_all_ic(1:nCluster,:)=actBSv;
%             actBSv_all_ic(nCluster+nM+1:nCluster_all_ic,nM+1:nB) = eye(nP);
% 
%             % Sj
%             Sj_rho_all_ic           =   zeros(nB, nCluster_all_ic);
%             Sj_rho_all_ic(:,1:nCluster) =   Sj_rho;
%             Sj_rho_all_ic(:,nCluster+1:nCluster_all_ic) = Sj_rho(:, BS_index);
% 

       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Get SINR from individual BSs
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [c_zf_BS, S_zf_BS, IN_zf_BS, SINR_zf_BS, c_mrt_BS, S_mrt_BS, IN_mrt_BS, SINR_mrt_BS] = SINR_BS_het(DistBU, nM, nP, L0, alpha_macro, alpha_pico, M, SU, Pt);
        fprintf('finish peak rate in cellular case for static users\n');


            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% save data
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf('end loop=%d, rho=%d\n', loop, step_rho);            
        save(strcat('data/','rate-',int2str(loop),'-rho-0',int2str(step_rho*100),'-nphot-',int2str(denseP_hot/denseP_reg),'-nuhot-',int2str(denseu_hot/denseu_reg),'-grid.mat'),'c_zf_rho','SU_rho','actBSv','nCluster','archv','narch','c_zf_orth','SU_rho_orth','actBSv_orth','archv_orth','c_zf_off','c_zf_BS');

        % =================================================
        save('data/loop.mat','loop');
        
end
            
