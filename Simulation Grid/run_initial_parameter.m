%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% generate multiple deployment realizations for simulations
% Last update 1/11/15 9:50pm

% save the deployments of BSs and users, as well as the cvx results
% users are static
% small BSs are uniformly dropped

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; clf; close all;

R_movev                     =   [100];%50 100];
step_rhov                   =   [1 0.75 0.5 0.25];%, 0.75, 0.5, 0.25];%[(16/4-1)/3];  


for lr = 1 : length(R_movev)
    R_move                  =   R_movev(lr);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% run data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    loop                =   0;
    save('data/loop.mat','loop');
    
    load(strcat('data/','Case1-Rmove-',int2str(R_move),'-rho-0',int2str(1*100),'.mat'));
    while loop < totloop
        
        load('data/loop.mat');
        loop=loop+1;
        if mod(loop, 4) ==0
            clc;
        end
        fprintf('start loop=%d\n',loop);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Network Deployment
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % =================================================
        % Macro BS location
        BS                      =   zeros(nM+nP, 2);
        % hexgonal model: number of layer nL
        BStp                    =   hextop(nL*2-1, nL*2-1);
        
        
        
        %% Deployment, BSs, users,
        BS      = zeros(nm,2);
        BS(2,:) = [0, sqrt(3)*rhex];
        BS(3,:) = [3/2*rhex, sqrt(3)/2*rhex]; 
        BS(4,:) = [1, -1].*BS(3,:);
        BS(5,:) = [0, -sqrt(3)*rhex];
        BS(6,:) = [-1, -1].*BS(3,:);
        BS(7,:) = [-1, 1].*BS(3,:);
        BS(8,:) = [0, 2*sqrt(3)*rhex];
        BS(9,:)= [3/2*rhex, sqrt(3)*3/2*rhex]; 
        BS(10,:)= [3*rhex, sqrt(3)*rhex]; 
        BS(11,:)= [3*rhex, 0]; 
        BS(12,:)= [1, -1].*BS(10,:);
        BS(13,:)= [1, -1].*BS(9,:);
        BS(14,:)= [-1,-1].*BS(8,:);
        BS(15,:)= [-1,-1].*BS(9,:);
        BS(16,:)= [-1,-1].*BS(10,:);
        BS(17,:)= [-1, 1].*BS(11,:);
        BS(18,:)= [-1, 1].*BS(10,:);
        BS(19,:)= [-1, 1].*BS(9,:);
        BS(20,:)= [3/2*rhex, sqrt(3)*5/2*rhex]; 
        BS(21,:)= [3*rhex, 2*sqrt(3)*rhex];
        BS(22,:)= [1, -1].*BS(20,:); 
        BS(23,:)= [1, -1].*BS(21,:); 
        BS(24,:)= [-1,-1].*BS(20,:); 
        BS(25,:)= [-1,-1].*BS(21,:); 
        BS(26,:)= [-1, 1].*BS(20,:); 
        BS(27,:)= [-1, 1].*BS(21,:); 
        BS(28,:)= [0, sqrt(3)*3*rhex]; 
        BS(29,:)= [0, -sqrt(3)*3*rhex]; 

        %BS(1:nM, :)             =   getCor(sqrt(nM), L);

        % =================================================
        % Pico BS locations and User locations per macrocell

        % grid model
        % [BS(nM+1:nM+nP,:), BSlocation] =   getCor(sqrt(nP),L); 
        % randomly drop
        %BS(nM+1:nM+nP,:)    	=   rand(nP,2)*L - L/2; 


        for lm = 1 : nM
            % location of pico BSs
            [BS((lm-1)*nPperM+nM+1:lm*nPperM+nM, :), ~]   =   PicoDeploy(BS(lm,:), L/sqrt(nM), nPreg/nM, nPhot/nM);

            % location of users
            [U((lm-1)*nuperM+1:lm*nuperM,:), ~]   	=   UEdeploy(BS(lm,:), L/sqrt(nM), denseu_reg, denseu_hot);
        end

        % =================================================
        % Users are moving
        Um                      =   zeros(nu, 2);
        theta                   =   rand(nu, 1) * (2 * pi);
        r                       =   R_move *sqrt(rand(nu, 1));
        % determine whether the users are moving or not
        Pmove_tp                =   rand(nu,1);
        move_index              =   Pmove_tp < P_move;
        Um_index                =   find(move_index==1);
        Um(:, 1)                =   U(:, 1) + r .* cos(theta) .* move_index;
        Um(:, 2)                =   U(:, 2) + r .* sin(theta) .* move_index;
        % wrap around the users move outside (L*L area)
        Utp                     =   find(Um(:,1)>L/2);
        Um(Utp,1)               =   -L/2 + Um(Utp,1)-L/2;
        Utp                     =   find(Um(:,1)<-L/2);
        Um(Utp,1)               =   L/2 + Um(Utp,1)+L/2;
        Utp                     =   find(Um(:,2)>L/2);
        Um(Utp,2)               =   -L/2 + Um(Utp,2)-L/2;
        Utp                     =   find(Um(:,2)<-L/2);
        Um(Utp,2)               =   L/2 + Um(Utp,2)+L/2;
        
        % Distance between UE and BS (wrap round)
        DistBU                  =   BS_UE_Distance_Calculation(BS, U, L);
        % Distance between UE and BS after moving
        DistBU_m                =   BS_UE_Distance_Calculation(BS, Um, L);

        % =================================================
        % different rho
        for ls =   1 : length(step_rhov) 
            step_rho            =   step_rhov(ls);
        
            %load(strcat('data/','Case2-gridPico-Rmove-',int2str(R_move),'-rho-0',int2str(step_rho*100),'.mat'));
            load(strcat('data/','Case1-Rmove-',int2str(R_move),'-rho-0',int2str(step_rho*100),'.mat'));

            
            % =================================================
            % MIMO parameters                       
            % user factor
            rho                     =   ones(nB, 1);
            for i = 2 : nB
                rho(i)              =   step_rho * i;
            end
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Get SINR from all possible clusters, with maximial cluster size NC_max
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
            % =================================================
            
            
            
            % get SINR, where each BS in cluster with size n can serve nS users
            [c_zf_rho, ~, ~, SINR_zf_rho, c_mrt_rho, S_mrt, IN_mrt, SINR_mrt, Sj_rho, actBSv, archvtp, nCluster] = SINR_het_rho_nS(NC_max, NC_min, DistBU, nM, nP, L0, alpha_macro, alpha_pico, M, Uj, Pt, rho);
            fprintf('finish peak rate for static users, nCluster=%d\n', nCluster);



            % =================================================
            % find the architecture-cluster indicator matrix
%             archv                   =   archvtp(NC_min:NC_max,:);
%             % number of architectures
%             narch                   =   size(archv, 1);


            % =================================================
            % Users are moving
            
            [c_zf_rho_m, ~, ~, SINR_zf_rho_m, c_mrt_rho_m, S_mrt_m, IN_mrt_m, SINR_mrt_m, ~, ~, ~,~] = SINR_het_rho_nS(NC_max, NC_min, DistBU_m, nM, nP, L0, alpha_macro, alpha_pico, M, Uj, Pt, rho);
            fprintf('finish peak rate for moving users, nCluster=%d\n', nCluster);

            % =================================================
            % find the index of macro BSs, pico BSs, and clusters with multiple BSs
%             macro_index             =   find(sum(actBSv(:,1:nM),2)==1 & sum(actBSv,2)==1);
%             pico_index              =   find(sum(actBSv(:,nM+1:nB),2)==1 & sum(actBSv,2)==1);
%             BS_index                =   [macro_index; pico_index];
%             cluster_index           =   find(sum(actBSv,2)>1);


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Get SINR when the macro BSs and pico are using orthogonal bands
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            % get SINR, where each BS in cluster with size n can serve nS users
            [c_zf_orthtp, ~, ~, SINR_zf_orth, c_mrt_orth, S_mrt_orth, IN_mrt_orth, SINR_mrt_orth, Sj_rho_orthtp, actBSv_orthtp, archvtp_orthtp, nCluster_orthtp] = SINR_het_rho_nS(NC_max, NC_min, DistBU(nM+1:nB,:), 0, nP, L0, alpha_macro, alpha_pico, M(nM+1:nB), Uj(nM+1:nB), Pt(nM+1:nB), rho);
            fprintf('finish peak rate for static users in orthogonal case, nCluster=%d\n', nCluster_orthtp+nM);
            macro_index             =   find(sum(actBSv(:,1:nM),2)==1 & sum(actBSv,2)==1);
            pico_index              =   find(sum(actBSv(:,nM+1:nB),2)==1 & sum(actBSv,2)==1);
            BS_index                =   [macro_index; pico_index];
            c_zf_orth               =   [c_zf_rho(macro_index,:);c_zf_orthtp];
            
            narch_orth              =   narch + 1;
            nCluster_orth           =   nCluster_orthtp + nM;
            Sj_rho_orth             =   zeros(nB, nCluster_orth);
            Sj_rho_orth(1:nM,1:nM) =   Sj_rho(1:nM, macro_index);
            Sj_rho_orth(nM+1:nB,nM+1:nCluster_orth) = Sj_rho_orthtp;
            
            actBSv_orth             =   zeros(nCluster_orth, nB);
            actBSv_orth(1:nM,1:nM)  =   actBSv(macro_index, 1:nM);
            actBSv_orth(nM+1:nCluster_orth,nM+1:nB) = actBSv_orthtp;
            
            archv_orth              =   zeros(narch_orth, nCluster_orthtp);
            archv_orth(1,1:nM)      =   archv(1, macro_index);
            archv_orth(2:narch_orth,nM+1:nCluster_orth) = archvtp_orthtp(NC_min:NC_max,:);

            [c_zf_orthtp, ~, ~, ~,~, ~, ~, ~, ~, ~, ~, ~] = SINR_het_rho_nS(NC_max, NC_min, DistBU_m(nM+1:nB,:), 0, nP, L0, alpha_macro, alpha_pico, M(nM+1:nB), Uj(nM+1:nB), Pt(nM+1:nB), rho);
            c_zf_orth_m             =   [c_zf_rho_m(macro_index,:);c_zf_orthtp];
            
            
            
           
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Get SINR when the macro BSs are off
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           
            % =================================================
            % spectral efficiency when macro BSs are off
            [c_zf_off, S_zf_off, IN_zf_off, SINR_zf_off, c_mrt_off, S_mrt_off, IN_mrt_off, SINR_mrt_off] = SINR_BS_het(DistBU, nM, nP, L0, alpha_macro, alpha_pico, M, Uj, Pt_ic);
            fprintf('finish peak rate in eICIC case for static users\n');

            % spectral efficiency in arch2+eICIC
            %c_zf_all_ic         	=   [c_zf_rho;c_zf_off];

            % =================================================
            % spectral efficiency of moving users when macro BSs are off 
            [c_zf_off_m, S_zf_off_m, IN_zf_off_m, SINR_zf_off_m, c_mrt_off_m, S_mrt_off_m, IN_mrt_off_m, SINR_mrt_off_m] = SINR_BS_het(DistBU_m, nM, nP, L0, alpha_macro, alpha_pico, M, Uj, Pt_ic);
            fprintf('finish peak rate in eICIC case for moving users\n');

            % spectral efficiency of moving users in arch2+eICIC
            %c_zf_all_ic_m           =   [c_zf_rho_m;c_zf_off_m];

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
            %% Get SINR from individual BSs in orth cases
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % get SINR with only macro BSs
            Pt_macro                =   Pt;
            Pt_macro(1+nM:nB)       =   0;
            [c_zf_macro, S_zf_macro, IN_zf_macro, SINR_zf_macro, c_mrt_macro, S_mrt_macro, IN_mrt_macro, SINR_mrt_macro] = SINR_BS_het(DistBU, nM, nP, L0, alpha_macro, alpha_pico, M, Uj, Pt_macro);
            fprintf('finish peak rate in orth case for static users\n');
            [c_zf_macro_m, S_zf_macro_m, IN_zf_macro_m, SINR_zf_macro_m, c_mrt_macro_m, S_mrt_macro_m, IN_mrt_macro_m, SINR_mrt_macro_m] = SINR_BS_het(DistBU_m, nM, nP, L0, alpha_macro, alpha_pico, M, Uj, Pt_macro);

            
            c_zf_orth_BS            =   [c_zf_macro(1:nM,:);c_zf_off(1+nM:nB,:)];
            c_zf_orth_BS_m          =   [c_zf_macro_m(1:nM,:);c_zf_off_m(1+nM:nB,:)];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Get SINR from individual BSs
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [c_zf_BS, S_zf_BS, IN_zf_BS, SINR_zf_BS, c_mrt_BS, S_mrt_BS, IN_mrt_BS, SINR_mrt_BS] = SINR_BS_het(DistBU, nM, nP, L0, alpha_macro, alpha_pico, M, Uj, Pt);
            fprintf('finish peak rate in cellular case for static users\n');

            % =================================================
            % Users are moving
            [c_zf_BS_m, S_zf_BS_m, IN_zf_BS_m, SINR_zf_BS_m, c_mrt_BS_m, S_mrt_BS_m, IN_mrt_BS_m, SINR_mrt_BS_m] = SINR_BS_het(DistBU_m, nM, nP, L0, alpha_macro, alpha_pico, M, Uj, Pt);
            fprintf('finish peak rate in cellular case for moving users\n');


            
            
            fprintf('end loop=%d, rho=%d\n',loop, step_rho);            
            save(strcat('data/',int2str(loop),'-rho-0',int2str(step_rho*100),'-deploy1.mat'),'BS','U','Um','c_zf_rho','Sj_rho','actBSv','nCluster','archv','narch','c_zf_rho_m','c_zf_orth','c_zf_orth_m','c_zf_orth_BS','c_zf_orth_BS_m','Sj_rho_orth','actBSv_orth','archv_orth','c_zf_off','c_zf_off_m','c_zf_BS','c_zf_BS_m');
            %save(strcat('data/',int2str(loop),'-rho-0',int2str(step_rho*100),'-deploy2.mat'),'BS','U','Um','c_zf_rho','Sj_rho','actBSv','nCluster','archv','narch','c_zf_rho_m','c_zf_orth','c_zf_orth_m','c_zf_orth_BS','c_zf_orth_BS_m','Sj_rho_orth','actBSv_orth','archv_orth','c_zf_off','c_zf_off_m','c_zf_BS','c_zf_BS_m');

            
            save('data/loop.mat','loop');
            end
    end
end
            
