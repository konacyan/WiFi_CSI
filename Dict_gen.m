function [PhiU, tempU, PhiA] = Dict_gen(Dicthe, Dictau, M, d, lambda, PU)
%% parameter
% Cortemplate:template of correlation matrix
% lags:length of correlation matrix
% fs:sampling frequency [63.9GHz ranging sampling clock]
% y = rcosdesign(beta,span,sps);
% beta: roll-off factor [0,1]
% dic_aoa: resolution of aoa estimation

%% Template
[I_1, ~]  = size(PU);

tempV = fix(ceil(Dictau)/I_1);
tempU = unique(ceil(Dictau) - tempV*I_1);

% PhiU: dictionary of the signal correlation vector
PhiU = PU(:,tempU);
% PhiV: dictionary of the delay 


% PhiA
A = @(theta)steering_vec(theta/180*pi,M,d,lambda);
% PhiA: dictionary of the array response
PhiA = A(Dicthe);

end