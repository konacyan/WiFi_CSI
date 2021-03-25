clc
clear 
% only have LOS path for test version

%% load data
addpath('tensorlab_2016-03-28');
filepath = 'test1_12.dat';
WLAN = wifi_config();
csi_data = gene_csi(filepath);
save('test1_12.mat','csi_data')

csi_length = length(csi_data);
for i_csi = 1 : csi_length
    csi = csi_data{i_csi};
end

% %----------test----------------- 
% csi = zeros(2,3,30);
% for i = 1:2
%     for j = 1:3
%         for k = 1:30
%             csi(i,j,k) = i + 0.1*j +0.001*k;
%             
%         end
%     end
% end
% %----------test----------------- 
%% generate simulation CSI57 data
% system setting
fc = 5.8e9; % carrier signal frequency [5.725-5.875] GHz
f_spacing = 312.5e3;% subcarrier spacing
c = 3e8; % speed of light
lambda = c/fc;  % wavelength(m)
d = lambda/2;
A = @(theta,M)steering_vec(theta/180*pi,M,d,lambda);

% multipath parameter
L = 3;% multhpath
theta = 32;% degree
dis = 7;% meter
sigma_tau = 3;% meter
ap_pos = [0,0];
agent_pos = [dis*cos(theta./180*pi),dis*sin(theta./180*pi)];
% add more constrains on the difference between ob1 and ob2?
ob1_pos = agent_pos + [-1*abs(sigma_tau*randn(1,1)), sigma_tau*randn(1,1)]; % obstacal1 postion
ob2_pos = agent_pos + [-1*abs(sigma_tau*randn(1,1)), sigma_tau*randn(1,1)]; % obstacal2 postion
Theta = [theta,180/pi*atan((ob1_pos(2)-ap_pos(2))/(ob1_pos(1)-ap_pos(1))),180/pi*atan((ob2_pos(2)-ap_pos(2))/(ob2_pos(1)-ap_pos(1)))];% -90~90 degree
Phi = -1.*[theta,180/pi*atan((ob1_pos(2)-agent_pos(2))/(ob1_pos(1)-agent_pos(1))),180/pi*atan((ob2_pos(2)-agent_pos(2))/(ob2_pos(1)-agent_pos(1)))]; % -90~90 degree
Dis = [dis,sqrt(sum((ap_pos - ob1_pos).^2)) + sqrt(sum((ob1_pos - agent_pos).^2)),sqrt(sum((ap_pos - ob2_pos).^2)) + sqrt(sum((ob2_pos - agent_pos).^2))];
Tau = Dis./c;

% sigma_theta = 10;%degree
% sigma_tau = 2;% meter
% Theta = repmat(theta,1,L)+sigma_theta*[0 randn(1,L-1)];% degree
% Dis = repmat(dis,1,L)+abs(sigma_tau*[0 randn(1,L-1)]); %meter
% Tau = Dis./c; %second

% generate csi data
N_t = 2;
N_r = 3;
N_c = 30; 
csi = zeros(N_t,N_r,N_c);
A_t = A(Phi,N_t);
A_r = A(Theta,N_r);
subcarrier_idx = [-28:2:-2,-1,1:2:27,28];% 20MHT mode, carriers index in IEEE 802.11n
f_carrier = fc + f_spacing.*subcarrier_idx;
Omega_tau = zeros(L,N_c);
for i_l = 1:1
    Omega_tau(i_l,:) = exp(-1j * 2 * pi .* f_carrier .* Tau(i_l));
    for i_t = 1:N_t
        for i_r = 1:N_r
            for i_c = 1:N_c
                csi(i_t,i_r,i_c) = A_t(i_t,i_l).*A_r(i_r,i_l).*Omega_tau(i_l,i_c);
            end
        end
    end
end


% add noise
SNR = 0; % dB
snr = 10^(SNR/10); % sigma_s^2
SINR = 10; % dB
% sinr = (snr./10.^( kron(SINR,ones(1,Num_test))/10 ))/L; % sigma_l^2 of each MPC
sinr = (snr./10.^( SINR/10 ))/L; % sigma_l^2 of each MPC




%% spatial smoothing
% %  original csi data([N_t, N_r, N_c])
% %  smoothing csi data([N_t * (1 + n_t_c),N_r * (1 + n_r_c), (N_c - n_t_c - n_r_c)])
% [N_t, N_r, N_c] = size(csi);
% n_t_c = 3; % number of virtual transmit array
% n_r_c = 3; % number of  virtual receiving array
% N_c_tilde = N_c - n_t_c - n_r_c;
% N_t_tilde = N_t * (1 + n_t_c); % original physical array plus virtual array
% N_r_tilde = N_r * (1 + n_r_c); % original physical array plus virtual array
% 
% csi_smooth = zeros(N_t_tilde, N_r_tilde, N_c_tilde);
% csi_idx = repmat(1:(n_r_c+1),1,N_t);
% for i_t = 1: N_t_tilde
%    
%     csi_slice = reshape(csi(ceil(i_t/(n_t_c+1)),:,:), N_r, N_c);
%     CSI_slice = zeros(N_r_tilde, N_c_tilde);
%     % smoothing
%     i_idx = 1;
%     for i_r = 1: N_r
%         for i_nrc = 1:(1 + n_r_c)
%             CSI_slice(i_idx,:) = csi_slice(i_r ,(i_nrc + csi_idx(i_t) -1 ):(N_c_tilde - 1 + i_nrc + csi_idx(i_t) -1));
%             i_idx = i_idx +1;
%         end      
%     end
%     csi_smooth(i_t,:,:) = reshape(CSI_slice,1,N_r_tilde,N_c_tilde);
% end

%% tensor decomposition
% how to estimate the exact number of multipath L in indoor environment?
% only have LOS path for test version
X = csi;
L = 1;
% X = csi_smooth;
X_decom = cpd(X,L); 
Tauhat = zeros(1,L);
Phihat = zeros(1,L);
Thetahat = zeros(1,L);
for i_l=1:L
    % AOD estimation
    Phihat(i_l) = MUSIC(X_decom{1,1}(:,i_l), N_t, d, lambda); 
    % AOA estimation
    Thetahat(i_l) = MUSIC(X_decom{1,2}(:,i_l), N_r, d, lambda); 
     
%     % delay estimation
%     [indu, indv] = T1_toa(X_decom{1,1}(:,i_l), X_decom{1,2}(:,i_l), PU, PV);
%     Tauhat(i_l) = (indv-1)*I_1 + indu-1;
%     iv = fix(ceil(Tau(i_l)*fs)/I_1);
%     if indv < iv
%         Tauhat(i_l) = inf;
%     end
end
