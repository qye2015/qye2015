%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% generate multiple deployment realizations for simulations
% Last update 1/27/15 1:50pm

% save the deployments of BSs and users:
% 1. macro BSs are deployed according to grid model
% 2. users are static
% 3. small BSs are uniformly dropped in both hotspots (3 picos) and regular area (1 pico)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; clf; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% network parameters initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of totloops
totloop                     =   100;
% =================================================
% network layout parameters
% inter-BS distance
L_interBS                	=   1000;
% areas to be simulated
L                           =   2000;
% ------------------------------------------
% density of macro BSs
denseM                      =   1 / L_interBS / L_interBS; 
% number of total macro BSs
nM                          =   round(denseM * L * L);
% ------------------------------------------
% density of pico BSs in regular area
denseP_reg                  =   4 / L_interBS / L_interBS;
% density of pico BSs in hotspots
denseP_hot                  =   3 * denseP_reg;
% number of pico BSs per regular area
nPreg                       =   round(denseP_reg * L_interBS / 2 * L_interBS / 2);
% number of pico BSs per hotspot area
nPhot                       =   round(denseP_hot * L_interBS / 2 * L_interBS / 2);
% number of total pico BSs
nP                          =   (2 * nPreg + 2 * nPhot) * nM;
% number of total BSs
nB                          =   nM + nP;
% ------------------------------------------
% density of users in regular area
denseu_reg                  =   15 * denseP_reg;
% density of users in hotspots
denseu_hot                  =   8 * denseu_reg;
% average user density 
denseU                      =   denseu_reg / 2 + denseu_hot / 2;
% number of users per regular area
nureg                       =   round(denseu_reg * L_interBS / 2 * L_interBS / 2);
% number of users per hotspot area
nuhot                       =   round(denseu_hot * L_interBS / 2 * L_interBS / 2);
% number of users 
nu                          =   (2 * nureg + 2 * nuhot) * nM;

% =================================================
% channel model
% number of macro BSs antennas
M_macro                     =   100;
% number of pico BSs antennas
M_pico                      =   40;
% ------------------------------------------
% number of users per macro BS
S_macro                     =   10;
% number of users per pico BS
S_pico                      =   4;
% ------------------------------------------
% maximial cluster size for the most general case
NC_max                  	=   4;
% minimal cluster size for the most general case
NC_min                      =   1;
% ------------------------------------------
% noise KTB (B=10MHz): -104dBm
noise                       =   1.38 * 10^(-23) * 290 * 10 * 10^(6);
% transmit power of macro BSs (46dBm) normalized by noise
P_macro                     =   10 .^ (4.6) / 1000 / noise;
% transmit power of pico BSs (35dBm) normalized by noise
P_pico                      =   10 .^ (3.0) / 1000 / noise;
% ------------------------------------------
% path loss exponent of macro BSs
alpha_macro                 =   3.5;
% path loss exponent of pico BSs
alpha_pico                  =   4;
% reference distance
L0                          =   40;
% ------------------------------------------
% biasing factor of BSs
bias                        =   ones(nB, 1);

% =================================================
% limit of the number of iterations
iter_max_dualgrad           =   5 * 10^3;
% limit of the number of iterations for vq
max_vq                      =   40000;
max_loop_ts                 =   5000;

% =================================================
% Macro BS locations
BS                       	=   getCor(sqrt(nM), L);
% initialize for user locations
U                           =   zeros(nu, 2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% run realizations 
%(small cells and users are randomly deployed at each hotspot and regular area)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
loop                        =   0;
save('data/loop.mat','loop');
% =================================================
while loop < totloop
    
    % =================================================
    load('data/loop.mat');
    loop                    =   loop + 1;
    if mod(loop, 5) ==0
        clc;
    end
    fprintf('start loop=%d\n',loop);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% deploy small BSs
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % deploy 3 picos per hotspots uniformly and 1 pico per regular area
    % deploy users 10*nureg per hotspots and nureg per regular area uniformly
    for lm = 1 : nM
        % =================================================
        % location of pico BSs
        picohot             =   rand(nPhot * 2, 2) * L_interBS .* [ones(nPhot, 2); -ones(nPhot, 2)];
        picoreg             =   [-L_interBS / 2, L_interBS / 2; L_interBS / 2, -L_interBS / 2];
        % pico index
        pico_index_tp    	=   nM + 1 + (lm-1) * (nPhot + nPreg) * 2 :  nM + lm * (nPhot + nPreg) * 2;  
        BS(pico_index_tp, :)=   [picohot; picoreg];
        BS(pico_index_tp, 1)=   BS(pico_index_tp, 1)+ BS(lm, 1);
        BS(pico_index_tp, 2)=   BS(pico_index_tp, 2)+ BS(lm, 2);
        
        % =================================================
        % location of users
        uehot               =   rand(nuhot * 2, 2) * L_interBS .* [ones(nuhot, 2); -ones(nuhot, 2)];
        uereg             	=   rand(nureg * 2, 2) * L_interBS .* [-ones(nureg, 1), ones(nureg, 1); ones(nureg, 1), -ones(nureg, 1)];
        % pico index
        ue_index_tp       	=   (lm-1) * (nuhot + nureg) * 2 + 1 : lm * (nuhot + nureg) * 2;  
        U(ue_index_tp, :) 	=   [uehot; uereg];
        U(ue_index_tp, 1)   =   U(ue_index_tp, 1) + BS(lm, 1);
        U(ue_index_tp, 2)   =   U(ue_index_tp, 2) + BS(lm, 2);
    end
    
    % =================================================
    % Distance between UE and BS (wrap round)
    DistBU              	=   BS_UE_Distance_Calculation(BS, U, L);
    
    % =================================================
    fprintf('end loop=%d', loop); 
    
    % =================================================
	save(strcat('data/',int2str(loop),'-nphot-',int2str(denseP_hot/denseP_reg),'-nuhot-',int2str(denseu_hot/denseu_reg),'-grid.mat'),'totloop','L_interBS','L','denseM','nM','denseP_reg','denseP_hot','nPreg','nPhot','nP','nB','denseu_reg','denseu_hot','denseU','nureg','nuhot','nu','M_macro','M_pico','S_macro','S_pico','NC_max','NC_min','noise','P_macro','P_pico','alpha_macro','alpha_pico','L0','bias','iter_max_dualgrad','max_vq','max_loop_ts','BS','U','DistBU');
    save('data/loop.mat','loop');
    
end
            
