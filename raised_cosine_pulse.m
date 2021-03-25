function y = raised_cosine_pulse(K, beta, span, sps, fs, tau)
% beta: roll-off factor [0,1]
% span: number of symbols [3]
% sps: samples per symbol= fb
% t: continuous time
% K: pulse duration  [2.4ns]
% t=k/fb where  k \in [-K/2:K/2]

%% signal generation
y = zeros(1,K);
% y = rcosdesign(beta,span,sps);
y(1+ceil(tau*fs):1+ceil(tau*fs)+span*sps) = rcosdesign(beta,span,sps);

% y = 1*sinc((k/fb-tau)/T).*cos(pi*beta*(k/fb-tau)/T)./(1 - (2*beta*(k/fb-tau)/T).^2);
% y( T/(2*beta)*fb + K/2 + tau) = pi/(4)*sinc(1/(2*beta));
% y( -T/(2*beta)*fb + K/2 + tau) = pi/(4)*sinc(1/(2*beta));
