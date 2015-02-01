% Function to get the max-biased-sinr association
% 
% ===================================================
% Inputs:
%   SINR:   (nB*nU matrix)      the SINR from each BS to each user
%   bias:   (nB*1 vector)   	the bias factor of each BS
%   c:      (nB*nU matrix)    	spectral efficiency from BS to user
% ===================================================
% Outputs:
%   x:      (nB*nU matrix)      association indictor
%   act_frac:(nB*nU matrix)   	activity fraction
%   rate:   (1*nU vecotr)       rate of each user
%   utility: (number)           utility function
%
% last updated: 6/26/14 2:56pm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x, act_frac, rate, utility] = max_bias_SINR(SINR, bias, c, Uj)

% number of BSs
nB                  =   length(SINR(:, 1));
% number of users
nU                  =   length(SINR(1, :));
% =================================================
% find association
SINR_bias               =   SINR .* repmat(bias, 1, nU);
x                       =   (SINR_bias == repmat(max(SINR_bias, [], 1), nB, 1));
% normalized factor per BSs
factor                  =   repmat(Uj ./ sum(x, 2), 1, nU);
act_frac                =   x .* min(factor, 1);

% rate
rate                    =   sum(c.*act_frac, 1);

% get utility
utility                 =   sum(log(rate));

end








