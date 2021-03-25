function [PU, PV] = Tem_gen(I_1, I_2, Cortemplate, lags)
%% parameter
% Cortemplate:(vector) template of the correlation matrix
% lags:length of correlation matrix
% fs:sampling frequency [63.9GHz ranging sampling clock]
% y = rcosdesign(beta,span,sps);
% beta: roll-off factor [0,1]


%% Template
% PU: basis of the signal correlation vector
PU = zeros(I_1);
PU(:,1) = [Cortemplate;zeros(I_1-lags,1)];
temp = PU(:,1);
for i = 2:(I_1-lags)
    PU(:,i) = circshift(temp,i-1,1);
end
for i = (I_1-lags+1):I_1
    PU(:,i) = [zeros(i-1,1); Cortemplate(1:I_1-i+1)];
end
% PV: basis of the delay
PV = diag(ones(I_2,1));

end