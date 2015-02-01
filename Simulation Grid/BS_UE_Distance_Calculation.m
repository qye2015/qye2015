% Function to get the wrap around distance between UEs and BSs
% 
% ===================================================
% Inputs:
%   BS:     (nB*2 matrix)       the BS locations
%   U:      (nu*2 matrix)   	the UE locations
%   Lsize: 	(number)            size of the network (Lsize * Lsize)
% ===================================================
% Outputs:
%   D:      (nB*nU matrix)      the distance between BS and UE
%
% last updated: 7/28/14 2:56pm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DistanceM  =   BS_UE_Distance_Calculation(BS, U, Lsize)

% number of users
nu                  =   length(U(:, 1));
% number of BSs
nB                  =   length(BS(:, 1));

% ===================================================
% creat a nB * nU matrix
BS_repx             =   repmat(BS(:, 1), 1, nu);
BS_repy             =   repmat(BS(:, 2), 1, nu);
U_repx              =   repmat(U(:, 1)', nB, 1);
U_repy              =   repmat(U(:, 2)', nB, 1);

% ===================================================
% the x distance between UE and BS
x_dist              =   abs(BS_repx - U_repx);
% the y distance between UE and BS
y_dist              =   abs(BS_repy - U_repy);

% ===================================================
% modulo operation due to wrap-around
x_dist              =   min(x_dist, Lsize - x_dist);
y_dist              =   min(y_dist, Lsize - y_dist);

% ===================================================
DistanceM           =   sqrt(x_dist .^ 2 + y_dist .^ 2);
end