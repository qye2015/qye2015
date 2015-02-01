% Function to get the distance between UEs and BSs
% the model of macro BSs is hexogonal 
% 
% ===================================================
% Inputs:
%   BS:     (nB*2 matrix)       the BS locations
%   U:      (nu*2 matrix)   	the UE locations
%   nm:     (number)            the # of macro BSs
%   R: 	(number)                inter-cell distance / 2
% ===================================================
% Outputs:
%   D:      (nB*nU matrix)      the distance between BS and UE
%  
% last updated: 1/21/15 4:56pm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function D = Hex_BS_UE_Distance_wraparound(BS, U, nm, R)

tp              =   hextop(10,10);

plot(tp(1,:),tp(2,:),'bo')




end