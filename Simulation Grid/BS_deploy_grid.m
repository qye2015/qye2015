% Function to deploy macro and pico BSs
% a center 7-hex-cell and wrap around
% 
% ===================================================
% Inputs:
%   BS:     (nB*2 matrix)       the BS locations
%   nm:     (number)            # of macro BSs
%   
%   L: 	(number)            # of layers, where each layer is consisted of a 7-hex-cell
% ===================================================
% Outputs:
%   D:      (nB*nU matrix)      the distance between BS and UE
%
% last updated: 1/26/15 2:56pm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DistanceM  =   BS_deploy_grid(BS, nm, L, )

% the hexgonal BS locations
BStp                =   getCor(sqrt(nP),L)

BS_index(5)         =   6;
BS_index(6)         =   1;





end