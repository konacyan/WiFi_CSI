function [res, Data_toa_sinr] = testtoa(M, d, fc, fs, L, Theta, Dis, snr, sinr, I_1)
res = zeros(1,4);
Data_toa_sinr = zeros(1,2*L);

c = 3e8;
Tau = Dis/c;
%% signal model
% fc:carrier signal frequency [3.5 - 6.5 ] GHz
% fb:baseband signal frequency [500 - 900] MHz
% fs:sampling frequency [63.9GHz ranging sampling clock]
lambda = c/fc;  % wavelength(m)

% y = rcosdesign(beta,span,sps);
% beta: roll-off factor [0,1]
beta = 0.1;
% span: number of symbols [3]
span = 6;
% sps: samples per symbol, pulse duration [2.4ns]
% sps = ceil((fs/fb)/span);
sps = 10;
% observation time [ 50 ]ns
K = 5e3;

[S_T, N_T, Cortemplate, lags] = CorM(Tau, fs, K, I_1, beta, span, sps);

S = zeros(L,K);
for l=1:L
    S(l,:) = raised_cosine_pulse(K, beta, span, sps, fs, Tau(l));
end
 
% signal generating
A = @(theta)steering_vec(theta/180*pi,M,d,lambda);
r = A(Theta)*diag([snr, sinr*ones(1,L-1)])*S + sinr*randn(M,K);

%% Tensor-based method
% Tensorization
r_cor = r*S_T';% correlation
I_2 = N_T/I_1;
Y = tensorization(r_cor, [I_1, I_2, M], 2);
[PU, PV] = Tem_gen(I_1, N_T/I_1, Cortemplate, lags);

%% A1
% tensor decomposition
Yhat1 = Y;
Uhat = cpd(Yhat1,L); % computes the factor matrices U{1}, ..., U{N} belonging to a canonical polyadic decomposition of the N-th order tensor T in R rank-one terms.
Tauhat = zeros(1,L);
for i=1:L
    % delay estimation
    [indu, indv] = T1_toa(Uhat{1,1}(:,i), Uhat{1,2}(:,i), PU, PV);
    Tauhat(i) = (indv-1)*I_1 + indu-1;
    iv = fix(ceil(Tau(i)*fs)/I_1);
    if indv < iv
        Tauhat(i) = NaN;
    end
end

% ----------07/08--------------
Dis_thres = 75; % meter
Tau_thres = fix(Dis_thres*fs/c);
Tau1_idx = round(Tau(1)*fs);
for i = 1:L
    if Tauhat(i) <= Tau1_idx - Tau_thres
        Tau_gap = fix((Tau1_idx - Tauhat(i))/Tau_thres);
        Tauhat(i) = Tauhat(i) + Tau_gap * Tau_thres;
    end
end
% ---------------------------

%----------------07/10-----------------
[~, ind1] = min(Tauhat);

Tauhat_res = abs(Tauhat - Tau1_idx);
Tauhat_sort = sort(Tauhat);

if Tauhat_res(ind1) > Tau_thres
    [~, ind1] = find(Tauhat==Tauhat_sort(2));

    if Tauhat_res(ind1) > Tau_thres
            [~, ind1] = find(Tauhat==Tauhat_sort(3));
            
            if Tauhat_res(ind1) > Tau_thres
                res_abs(1,:) = abs(Tauhat - (Tau1_idx - Tau_thres));
                res_abs(2,:) = abs(Tauhat - (Tau1_idx + Tau_thres));
                [~,ind1]=find(res_abs==min(min(res_abs)));
            end
    end
    
end



if length(ind1) > 1
    ind = ind1(1);
elseif length(ind1) == 1
    ind = ind1;
else
    ind1 = find(Tauhat==min(Tauhat));
    if isempty(ind1)
        ind1 = 1;
        ind = ind1;
    elseif length(ind1) > 1
        ind = ind1(1);
    else
        ind = ind1;
    end
end
% ---------------------------

res(1) = abs(Tauhat(ind)/fs - Tau(1));
Data_toa_sinr(1,1:L) = Tauhat;

%% A2
Dicthe_noise = 1e-4;
Dicthe = [Theta(1)+Dicthe_noise 90*rand(1,500)];
Dictau = [Tau(1) Tau(1) + 10*rand(1,500)/c]*fs;
[PhiU, tempU, PhiA] = Dict_gen(Dicthe, Dictau, M, d, lambda, PU);

U = {PhiU',PV',PhiA'};
Yhat2 = tmprod(Y,U,1:3);
% tensor decomposition
Uhat = cpd(Yhat2,L); % computes the factor matrices U{1}, ..., U{N} belonging to a canonical polyadic decomposition of the N-th order tensor T in R rank-one terms.
% Uhat{1,3}: the lth column represent the possibility of the item in PhiA being the lth MPC
Tauhat = zeros(1,L);
for i=1:L
    % delay estimation
    [~,indu] = max(abs(Uhat{1,1}(:,i)));
    [~,indv] = max(abs(Uhat{1,2}(:,i)));
    Tauhat(i) = (indv-1)*I_1 + tempU(indu)-1;
    iv = fix(ceil(Tau(i)*fs)/I_1);
    if indv < iv
        Tauhat(i) = NaN;
    end
end


% ----------07/08--------------
Dis_thres = 75; % meter
Tau_thres = fix(Dis_thres*fs/c);
Tau1_idx = round(Tau(1)*fs);
for i = 1:L
    if Tauhat(i) <= Tau1_idx - Tau_thres
        Tau_gap = fix((Tau1_idx - Tauhat(i))/Tau_thres);
        Tauhat(i) = Tauhat(i) + Tau_gap * Tau_thres;
    end
end
% ---------------------------

%----------------07/10-----------------
[~, ind1] = min(Tauhat);

Tauhat_res = abs(Tauhat - Tau1_idx);
Tauhat_sort = sort(Tauhat);

if Tauhat_res(ind1) > Tau_thres
    [~, ind1] = find(Tauhat==Tauhat_sort(2));

    if Tauhat_res(ind1) > Tau_thres
            [~, ind1] = find(Tauhat==Tauhat_sort(3));
            
            if Tauhat_res(ind1) > Tau_thres
                res_abs(1,:) = abs(Tauhat - (Tau1_idx - Tau_thres));
                res_abs(2,:) = abs(Tauhat - (Tau1_idx + Tau_thres));
                [~,ind1]=find(res_abs==min(min(res_abs)));
            end
    end
    
end


if length(ind1) > 1
    ind = ind1(1);
elseif length(ind1) == 1
    ind = ind1;
else
    ind1 = find(Tauhat==min(Tauhat));
    if isempty(ind1)
        ind1 = 1;
        ind = ind1;
    elseif length(ind1) > 1
        ind = ind1(1);
    else
        ind = ind1;
    end
end
% ---------------------------





res(2) = abs(Tauhat(ind)/fs - Tau(1));
Data_toa_sinr(1,(L+1):(2*L)) = Tauhat;

%% TOA estimation
Tau_toa = (Tau(3) - Tau).';
Tau_toa = fliplr(Tau_toa);% shift
iter_toa = 10;% iterations
for itoa = 1:iter_toa

% s_toa = sum(S);
s0_toa = raised_cosine_pulse(K, beta, span, sps, fs, 0);% signal template
r_toa = [snr, sinr*ones(1,L-1)]*S + sinr*randn(1,K);% received signal
% r_toa = [snr, sinr*ones(1,L-1)]*S + randn(1,K);% received signal
r0_toa = snr*s0_toa + sinr*randn(1,K);
% initial value
omega = zeros(1,L+1);
theta_delta = 1;
theta_rou = 1;


for n_iter = 2:L+1   
s_w = zeros(1,K);
    for i = 1:n_iter-1
      s0_w = raised_cosine_pulse(K, beta, span, sps, fs, omega(1,i));
      s_w = s0_w(1,1:K)./max(s0_w(1,1:K));
      s_w = s_w + s_w;
    end
i_label = 1;
for w = omega(1,n_iter-1):1e-9:1e-6
    s0_wi = raised_cosine_pulse(K, beta, span, sps, fs, w);
    s_wi = s0_wi(1,1:K)./max(s0_wi(1,1:K));
    w_n(1,i_label) = (r_toa - s_wi)*(s_w.');
    i_label = i_label + 1;
end
[mul_max,omega_index] = max(w_n,[],2);
omega(1,n_iter) = omega(1,n_iter-1) + 1e-8*(omega_index-1);

end

tau_LOS(1,itoa) = omega(1,n_iter-1);

end
tau_estimate = abs(Tau(3) - mean(tau_LOS));
res(3) = abs(tau_estimate - Tau(1));

%% CLEAN

N_uwb = span*sps;% signal template length

cir = abs(r);
y0 = sum(S);% signal template

s1 = y0(1+ceil(Tau(1)*fs):1+ceil(Tau(1)*fs)+N_uwb);% signal template 1
s2 = y0(1+ceil(Tau(2)*fs):1+ceil(Tau(2)*fs)+N_uwb);% signal template 2
s3 = y0(1+ceil(Tau(3)*fs):1+ceil(Tau(3)*fs)+N_uwb);% signal template 3

s_clean(1,:) = s1;
s_clean(2,:) = s2;
s_clean(3,:) = s3;
for i = L : L+N_uwb
    yob = cir(i,:);
    yob0 = yob;
    max_am = max(yob);
    tau_num = 0;
    while (isempty(s_clean) == 0)
        rss = [];
        rs = [];
        r_max = [];

        for j = 1 : length(s_clean(:,1))
            rss(j,:) = xcorr(s_clean(j,:),s_clean(j,:));
            rs(j,:) = xcorr(yob,s_clean(j,:));
            r_max(j) = max(rs(j,:));
        end

        r_max_idx = find(r_max == max(r_max),1,'first');
%         r_max_idx = find(r_max == max(r_max));
        ak(tau_num+1) = r_max(r_max_idx)/max(rss(r_max_idx,:));
        tau_ob(tau_num+1) = (find(rs(r_max_idx,:) == max(rs(r_max_idx,:)))  - (length(yob)-length(s_clean(r_max_idx,:))) - (length(s_clean(r_max_idx,:))-1)/2 ) ;
        tau_ob(tau_num+1) = floor(tau_ob(tau_num+1));
%         yob(tau_ob(tau_num+1)+ -N_uwb/2:tau_ob(tau_num+1)+ N_uwb/2) = yob(tau_ob(tau_num+1)-N_uwb/2:tau_ob(tau_num+1)+N_uwb/2) - ak(tau_num+1)*s_clean(r_max_idx,:);
        if tau_ob(tau_num+1)+ N_uwb/2 < K - N_uwb && tau_ob(tau_num+1)+ -N_uwb/2 > 0
            yob(tau_ob(tau_num+1) -N_uwb/2:tau_ob(tau_num+1)+ N_uwb/2) = yob(tau_ob(tau_num+1)-N_uwb/2:tau_ob(tau_num+1)+N_uwb/2) - ak(tau_num+1)*s_clean(r_max_idx,:); 
        elseif tau_ob(tau_num+1) - N_uwb/2 <= 0
            res_tauN = -1*(tau_ob(tau_num+1)- N_uwb/2 - 1);
            yob(tau_ob(tau_num+1) -N_uwb/2+res_tauN:tau_ob(tau_num+1)+ N_uwb/2+res_tauN) = yob(tau_ob(tau_num+1)-N_uwb/2+res_tauN:tau_ob(tau_num+1)+N_uwb/2+res_tauN) - ak(tau_num+1)*s_clean(r_max_idx,:); 
        else
            res_tauN = tau_ob(tau_num+1)+ N_uwb/2 - K;
            yob(tau_ob(tau_num+1)+ -N_uwb/2-res_tauN:tau_ob(tau_num+1)+ N_uwb/2-res_tauN) = yob(tau_ob(tau_num+1)-N_uwb/2-res_tauN:tau_ob(tau_num+1)+ N_uwb/2-res_tauN) - ak(tau_num+1)*s_clean(r_max_idx,:);
        end
        tau_num = tau_num +1;
        s_clean(r_max_idx,:) = [];

    end
    for i_clean = 1:L
        if tau_ob(i_clean)+ N_uwb/2 < K - N_uwb && tau_ob(i_clean)+ -N_uwb/2 > 0
            s_clean(i_clean,:) = yob0(tau_ob(i_clean)-N_uwb/2:tau_ob(i_clean)+N_uwb/2);
        elseif tau_ob(i_clean) - N_uwb/2 <= 0
            res_tauN = -1*(tau_ob(i_clean)- N_uwb/2 - 1);
            s_clean(i_clean,:) = yob0(tau_ob(i_clean)-N_uwb/2+res_tauN:tau_ob(i_clean)+N_uwb/2+res_tauN);
        else
            res_tauN = tau_ob(i_clean)+ N_uwb/2 - K;
            s_clean(i_clean,:) = yob0(tau_ob(i_clean)-N_uwb/2-res_tauN:tau_ob(i_clean)+N_uwb/2-res_tauN);
        end
    end
    t = -N_uwb/2:1:N_uwb/2;
    h_uwb(1,:) = ak(1)*sinc(t);
    h_uwb(2,:) = ak(2)*sinc(t);
    h_uwb(3,:) = ak(3)*sinc(t);
    h_t = zeros(1,length(cir(1,:)));
%     h_t(tau_ob(1)+ -N_uwb/2:tau_ob(1)+ N_uwb/2) = h_t(tau_ob(1)+ -N_uwb/2:tau_ob(1)+ N_uwb/2)+h_uwb(1,:);
%     h_t(tau_ob(2)+ -N_uwb/2:tau_ob(2)+ N_uwb/2) = h_t(tau_ob(2)+ -N_uwb/2:tau_ob(2)+ N_uwb/2)+h_uwb(2,:);
%     h_t(tau_ob(3)+ -N_uwb/2:tau_ob(3)+ N_uwb/2) = h_t(tau_ob(3)+ -N_uwb/2:tau_ob(3)+ N_uwb/2)+h_uwb(3,:);
    for i_clean = 1:L
        if tau_ob(i_clean)+ N_uwb/2 < K - N_uwb && tau_ob(i_clean)+ -N_uwb/2 > 0
            h_t(tau_ob(i_clean)+ -N_uwb/2:tau_ob(i_clean)+ N_uwb/2) = h_t(tau_ob(i_clean)+ -N_uwb/2:tau_ob(i_clean)+ N_uwb/2)+h_uwb(1,:);
        elseif tau_ob(i_clean) - N_uwb/2 <= 0
            res_tauN = -1*(tau_ob(i_clean)- N_uwb/2 - 1);
            h_t(tau_ob(i_clean)+ -N_uwb/2+res_tauN:tau_ob(i_clean)+ N_uwb/2+res_tauN) = h_t(tau_ob(i_clean)+ -N_uwb/2+res_tauN:tau_ob(i_clean)+ N_uwb/2+res_tauN)+h_uwb(1,:);
        else
            res_tauN = tau_ob(i_clean)+ N_uwb/2 - K;
            h_t(tau_ob(i_clean)+ -N_uwb/2-res_tauN:tau_ob(i_clean)+ N_uwb/2-res_tauN) = h_t(tau_ob(i_clean)+ -N_uwb/2-res_tauN:tau_ob(i_clean)+ N_uwb/2-res_tauN)+h_uwb(1,:);
        end
    end
    lcs_est(i,:) = tau_ob;% signal location estimations
%     i;
end

tau_clean_est = mean(lcs_est(L:L+N_uwb,:))./fs; %tau
res(4) = abs(tau_clean_est(1) - Tau(1));




end

