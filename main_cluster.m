clc
clear 
close all
addpath('tensorlab_2016-03-28');
addpath('linux-80211n-csitool-supplementary');
%% code parameter
simulation_flag = 0;
plot_flag = 1;
phaseCalibration_flag = 1;
bandwidth_flag = 0; %if bandwith flag == 1, HT40MHz mode; else HT20MHz mode
antenna_flag = 0; % if antenna_flag == 1, N_tx = 3; else N_tx = 2;

%% system setting
% fc = 2.472e9; % Channel 13
fc = 5.2e9; 
f_gap = 312.5e3;% subcarrier spacing
c = 3e8; % speed of light
lambda = c/fc;  % wavelength(m)
d = 5.35e-2./2;
A = @(theta,M)steering_vec(theta/180*pi,M,d,lambda);
if (bandwidth_flag)
    subcarrier_idx = -58:4:58;% 40MHT mode, carriers index in IEEE 802.11n
    subcarrier_spacing = 4 * f_gap;
else
    subcarrier_idx = [-28:2:-2,-1,1:2:27,28];% 20MHT mode, carriers index in IEEE 802.11n
    Subcarrier_idx = -30:2:28;
    subcarrier_spacing = 2 * f_gap;
end

f_carrier = fc + f_gap.*subcarrier_idx;

%% load data
filename = '0316log12_2';
filepath = ['Data/0316_AP/',filename,'.dat'];
WLAN = wifi_config();
[csi_data,csi_trace] = gene_csi(filepath);
% save([filename,'.mat'], 'csi_trace');

% Groundtruth
groundtruth_x = str2num(filename(end-3));
groundtruth_y = str2num(filename(end-2));
groundtruth =  180/pi*atan(groundtruth_y/groundtruth_x)


Num_Pack = 10; % package number
for i_csi = 1:Num_Pack
    csi_size =  size(csi_data{i_csi});
    csi = csi_data{i_csi};
    csi_detail = csi_trace{i_csi};
       
% sort csi
[N_t, N_r, N_c] = size(csi);
if csi_detail.perm == [1 3 2]
    csi_sort = zeros(N_t,N_r,N_c);
    csi_sort(:,1,:)=csi(:,1,:);
    csi_sort(:,3,:)=csi(:,2,:);
    csi_sort(:,2,:)=csi(:,3,:);
    csi = csi_sort;
end

%% plot figure
if plot_flag == 1
% plot SNR data
    csi_SNR_rxA = zeros(N_r,N_c);
    csi_SNR_rxB = zeros(N_r,N_c);
%     csi_SNR_rxC = zeros(N_r,N_c);
    for j=1:N_r   
        for k=1:N_c    
            csi_SNR_rxA(j,k)=csi(1,j,k);
            csi_SNR_rxB(j,k)=csi(2,j,k);
%             csi_SNR_rxC(j,k)=csi(3,j,k);
        end
    end
%     figure;
%     subplot(1,3,1)
%     plot(db(abs(csi_SNR_rxA.')))
%     legend('RX Antenna A', 'RX Antenna B', 'RX Antenna C');
%     xlabel('Subcarrier index');
%     ylabel('SNR [dB]')
%     subplot(1,3,2)
%     plot(db(abs(csi_SNR_rxB.')))
%     legend('RX Antenna A', 'RX Antenna B', 'RX Antenna C');
%     xlabel('Subcarrier index');
%     ylabel('SNR [dB]')
%     subplot(1,3,3)
%     plot(db(abs(csi_SNR_rxC.')))
%     legend('RX Antenna A', 'RX Antenna B', 'RX Antenna C');
%     xlabel('Subcarrier index');
%     ylabel('SNR [dB]')

% plot wrapped phase 
    N_pack = 1;% first package
    for i_t = 1:N_t
        for i_r = 1:N_r
            csi_pack = csi_data{N_pack};
            csi_subcarriers(:,((i_t - 1)*N_r + i_r)) = squeeze(csi_pack(i_t,i_r,:)).';
        end
    end
    csi_phase=angle(csi_subcarriers);
    csi_unwrapphase = unwrap(csi_phase,pi/2);
%     figure;
%     plot(180*csi_phase/pi,'.','Markersize',25);
%     xlabel('Subcarrier Index');
%     ylabel('Phase [deg]')

% plot unwrapped phase
%     figure;
%     for i_s = 1:N_r*N_t
%             subplot(N_t, N_r,i_s)
%             plot(csi_unwrapphase(:,i_s));
%             xlabel('Subcarrier Index');
%             ylabel('Unwrap Phase');
%             tx_num = fix((i_s-1)/N_r)+1;
%             rx_num = mod((i_s-1),N_r)+1;
%             title(['TX:',num2str(tx_num) '  RX:',num2str(rx_num)])
%     end

 
    % plot multipackage unwrapped phase 
%     figure;
%     N_pack = 20;% first package
% 
%     for i_n = 1:N_pack
%         
%         csi_multipack = csi_data{i_n};
%         csi_multisub(:,i_n) = squeeze(csi_multipack(1,1,:)).';
% 
%     end
%     csi_multiphase=angle(csi_multisub);
%     csi_multiunwrapphase = unwrap(csi_multiphase,pi/2);
%     plot(csi_multiunwrapphase);
%     title(filepath);
%     saveas(gcf,'Unwrap Phase','png');
    
    % plot the residual
%     figure;
%      for i_s = 1:N_r*N_t
%             subplot(N_t, N_r,i_s)
%             
%             tx_num = fix((i_s-1)/N_r)+1;
%             rx_num = mod((i_s-1),N_r)+1;
%             
%             plot(csi_unwrapphase(:,i_s) - csi_unwrapphase(:,((tx_num-1) * N_r + 1)) );
%             xlabel('Subcarrier Index');
%             ylabel('Unwrap Phase');
%             title(['TX:',num2str(tx_num) '  RX:',num2str(rx_num)])
%     end
    
    
end
%% phase interpolation
if phaseCalibration_flag == 1 
    
%     sum_coefficient = polyfit(subcarrier_idx, mean(csi_unwrapphase.',1),1);
%     sum_polyval = polyval(sum_coefficient,subcarrier_idx);
%     figure;    
%     plot(subcarrier_idx, mean(csi_unwrapphase.',1),'o',subcarrier_idx, sum_polyval,'-')
%     title('Sum poly');

% Raw Data
%     figure;
    coefficient = zeros(N_t*N_r,2);
    Coefficient = zeros(N_t*N_r,2);
    csi_poly = zeros(N_t*N_r,N_c);
    for i_t = 1: N_t*N_r

        coefficient(i_t,:) = polyfit(subcarrier_idx,csi_unwrapphase(:,i_t).',1);
        csi_poly(i_t,:) = polyval(coefficient(i_t,:),subcarrier_idx);
%         subplot(N_t,N_r,i_t)
%         plot(subcarrier_idx,csi_unwrapphase(:,i_t),'o',subcarrier_idx,csi_poly(i_t,:),'-');

    end
    
% Remove Singular Data   
    coefficient_k = coefficient(:,1);
    coefficient_ave = mean(coefficient_k);
    coefficient_std = std(coefficient_k);
%     coefficient_thres = 3*coefficient_std;
    coefficient_thres = 0.01;
    for i_co = 1:length(coefficient_k)
        if abs(coefficient_ave - coefficient_k(i_co) > coefficient_thres)
            coefficient_k(i_co) = 0;
        end
    end
    coefficient_k(find(coefficient_k == 0)) = [];
    coefficient_Ave = mean(coefficient_k);
    Coefficient(:,1) = coefficient_Ave;
    
% Fit Curve   
    ft = fittype(@(b,x) coefficient_Ave*x+b );
    for i_ft = 1:N_t*N_r
        
        fit_result = fit(subcarrier_idx.',csi_unwrapphase(:,i_ft),ft,'start',0);
        Coefficient(i_ft,2) = fit_result(0);
    end
end    

% Rebuild CSI Data
%     figure;
    if bandwidth_flag == 0 % subcarrier spacing is nonuniform

        
        for i_t = 1: N_t*N_r

            csi_Poly(i_t,:) = polyval(Coefficient(i_t,:),Subcarrier_idx);
%             subplot(N_t,N_r,i_t)
%             plot(subcarrier_idx,csi_unwrapphase(:,i_t),'o',Subcarrier_idx, csi_Poly(i_t,:),'-');

        end    
    end   
    
    csi_modified_phase = csi_unwrapphase - csi_Poly.';
    csi_STO = -1*coefficient(:,1)./(2*pi*subcarrier_spacing);% can be further used to get the multipath number(SpotFi, page 8)
%     figure;
%     for i_t = 1: N_t*N_r
%         subplot(N_t,N_r,i_t)
%         plot(Subcarrier_idx,csi_modified_phase(:,i_t));
%     end   
    
    
    csi_amp = abs(csi_subcarriers);
    csi_modified = csi_amp.*exp(1i.*(csi_Poly.'));
    csi_Modified = zeros(N_t, N_r, N_c);
    
    for i_t = 1: N_t
        for i_r = 1:N_r
            for i_c = 1:N_c
                csi_Modified(i_t,i_r,i_c) = csi_modified(i_c,(i_t-1)*N_r+i_r);
            end
        end
    end
    csi_unModified = csi;
    csi = csi_Modified;
    
%% phase calibration
% if phaseCalibration_flag == 1
%     figure;
%     coefficient = zeros(N_t*N_r,2);
%     csi_poly = zeros(N_t*N_r,N_c);
%     for i_t = 1: N_t*N_r
% 
%         coefficient(i_t,:) = polyfit(1:N_c,csi_unwrapphase(:,i_t).',1);
%         csi_poly(i_t,:) = polyval(coefficient(i_t,:),1:N_c);
%         subplot(N_t,N_r,i_t)
%         plot(1:N_c,csi_unwrapphase(:,i_t),'o',1:N_c,csi_poly(i_t,:),'-');
% 
%     end
% 
%     csi_modified_phase = csi_unwrapphase - csi_poly.';
%     csi_STO = -1*coefficient(:,1)./(2*pi*subcarrier_spacing);% can be further used to get the multipath number(SpotFi, page 8)
%     figure;
%     for i_t = 1: N_t*N_r
%         subplot(N_t,N_r,i_t)
%         plot(1:N_c,csi_modified_phase(:,i_t));
%     end
% 
% 
%     csi_amp = abs(csi_subcarriers);
%     csi_modified = csi_amp.*(cos(csi_modified_phase) + 1i*sin(csi_modified_phase));
%     csi_Modified = zeros(N_t, N_r, N_c);
%     for i_t = 1: N_t
%         for i_r = 1:N_r
%             for i_c = 1:N_c
%                 csi_Modified(i_t,i_r,i_c) = csi_modified(i_c,(i_t-1)*N_t+i_r);
%             end
%         end
%     end
%     csi_unModified = csi;
%     csi = csi_Modified;
% end

%% generate simulation CSI data 
if simulation_flag == 1
% multipath parameter
    L = 3;% multhpath
    theta = 32;% degree
    dis = 7;% meter
    sigma_tau = 2;% meter
    ap_pos = [0,0];
    agent_pos = [dis*cos(theta./180*pi),dis*sin(theta./180*pi)];
% add more constrains on the difference between ob1 and ob2?
    ob1_pos = agent_pos + [-1*abs(sigma_tau*randn(1,1)), sigma_tau*randn(1,1)]; % obstacal1 postion
    ob2_pos = agent_pos + [-1*abs(sigma_tau*randn(1,1)), sigma_tau*randn(1,1)]; % obstacal2 postion
    Theta = [theta,180/pi*atan((ob1_pos(2)-ap_pos(2))/(ob1_pos(1)-ap_pos(1))),180/pi*atan((ob2_pos(2)-ap_pos(2))/(ob2_pos(1)-ap_pos(1)))];% -90~90 degree
    Phi = -1.*[theta,180/pi*atan((ob1_pos(2)-agent_pos(2))/(ob1_pos(1)-agent_pos(1))),180/pi*atan((ob2_pos(2)-agent_pos(2))/(ob2_pos(1)-agent_pos(1)))]; % -90~90 degree
    Dis = [dis,sqrt(sum((ap_pos - ob1_pos).^2)) + sqrt(sum((ob1_pos - agent_pos).^2)),sqrt(sum((ap_pos - ob2_pos).^2)) + sqrt(sum((ob2_pos - agent_pos).^2))];
    Tau = Dis./c;

    Dis_sort = sort(Dis);
    Dis_idx = zeros(1,length(Dis));
    for i_dis = 1:length(Dis)
       Dis_idx(1,i_dis) = find(  Dis == Dis_sort(i_dis));
    end
% Wi-Fi parameter
    N_t = 2;
    N_r = 3;
    N_c = 30; 
    csi = zeros(N_t,N_r,N_c);
    A_t = A(Phi,N_t);
    A_r = A(Theta,N_r);


% add noise
    SNR = 0; % dB
    snr = 10^(SNR/10); % sigma_s^2
    SINR = 0; % dB
% sinr = (snr./10.^( kron(SINR,ones(1,Num_test))/10 ))/L; % sigma_l^2 of each MPC
    sinr = (snr./10.^( SINR/10 ))/L; % sigma_l^2 of each MPC
    amp = [snr, sinr*ones(1,L-1)]; % attenuation
% amp = [ones(1,L)]; % attenuation
    Omega_tau = zeros(L,N_c);

%generate CSI data
    for i_l = 1:L
        Omega_tau(i_l,:) = exp(-1j * 2 * pi .* f_carrier .* Tau(i_l));
        for i_t = 1:N_t
            for i_r = 1:N_r
                for i_c = 1:N_c
    %                  csi(i_t,i_r,i_c) =  csi(i_t,i_r,i_c) + A_t(i_t,i_l).*A_r(i_r,i_l).*Omega_tau(i_l,i_c);
                    csi(i_t,i_r,i_c) =  csi(i_t,i_r,i_c) + amp(i_l).*A_t(i_t,i_l).*A_r(i_r,i_l).*Omega_tau(i_l,i_c);
                end
            end
        end
    end
% csi = csi + 0.1*randn(N_t,N_r,N_c);
end
%% spatial smoothing
%  original csi data([N_t, N_r, N_c])
%  smoothing csi data([N_t * (1 + n_t_c),N_r * (1 + n_r_c), (N_c - n_t_c - n_r_c)])
[N_t, N_r, N_c] = size(csi);
n_t_c = 3; % number of virtual transmit array
n_r_c = 3; % number of  virtual receiving array
N_c_tilde = N_c - n_t_c - n_r_c;
N_t_tilde = N_t * (1 + n_t_c); % original physical array plus virtual array
N_r_tilde = N_r * (1 + n_r_c); % original physical array plus virtual array

csi_smooth = zeros(N_t_tilde, N_r_tilde, N_c_tilde);
csi_idx = repmat(1:(n_r_c+1),1,N_t);
for i_t = 1: N_t_tilde
   
    csi_slice = reshape(csi(ceil(i_t/(n_t_c+1)),:,:), N_r, N_c);
    CSI_slice = zeros(N_r_tilde, N_c_tilde);
    % smoothing
    i_idx = 1;
    for i_r = 1: N_r
        for i_nrc = 1:(1 + n_r_c)
            CSI_slice(i_idx,:) = csi_slice(i_r ,(i_nrc + csi_idx(i_t) -1 ):(N_c_tilde - 1 + i_nrc + csi_idx(i_t) -1));
            i_idx = i_idx +1;
        end      
    end
    csi_smooth(i_t,:,:) = reshape(CSI_slice,1,N_r_tilde,N_c_tilde);
end

%% SpotFi
L = 3;

% % sample CSI trace is a 90x1 vector where first 30 elements correspond to subcarriers for first rx antenna, second 30 correspond to CSI from next 30 subcarriers and so on.
% % replace sample_csi_trace with CSI from Intel 5300 converted to 90x1 vector
% sample_csi_traceTmp = load('sample_csi_trace');
% sample_csi_trace = sample_csi_traceTmp.sample_csi_trace;

% sample_csi_trace = reshape(csi(1,:,:),N_r,N_c).'; %这里问题很大
sample_csi_trace = csi; 

% N_c = length(subcarrier_idx); % number of subcarriers
% subCarrSize = 128;  % total number fo
T = N_t; % number of transmitter antennas

% MUSIC algorithm requires estimating MUSIC spectrum in a grid. 
% paramRange captures parameters for this grid
% For the following example, MUSIC spectrum is caluclated for 101 ToF (Time of flight) values spaced equally between -25 ns and 25 ns.
% MUSIC spectrum is calculated for for 101 AoA (Angle of Arrival) values between -90 and 90 degrees.
paramRange = struct;
paramRange.GridPts = [101 101 1]; % number of grid points in the format [number of grid points for ToF (Time of flight), number of grid points for angle of arrival (AoA), 1]
paramRange.delayRange = [-50 50]*1e-9; % lowest and highest values to consider for ToF grid. 
paramRange.angleRange = 90*[-1 1]; % lowest and values to consider for AoA grid.
do_second_iter = 0;
% paramRange.seconditerGridPts = [1 51 21 21];
paramRange.K = floor(N_r/2)+1; % parameter related to smoothing.  
paramRange.L = floor(N_c/2); % parameter related to smoothing.  
paramRange.T = floor(N_t/2)+1;
paramRange.deltaRange = [0 0]; 

maxRapIters = Inf;
useNoise = 0;
paramRange.generateAtot = 2;

% ToF sanitization code (Algorithm 1 in SpotFi paper)
% csi_plot = reshape(sample_csi_trace, N_c, N_r);

for i_t = 1:N_t
    csi_sanitization = squeeze(sample_csi_trace(i_t,:,:));
    csi_sanitization = reshape(csi_sanitization, N_c, N_r);
    [PhsSlope, PhsCons] = removePhsSlope(csi_sanitization,N_r,subcarrier_idx,N_c);
    ToMult = exp(1i* (-PhsSlope*repmat(subcarrier_idx(:),1,N_r) - PhsCons*ones(N_c,N_r) ));
    csi_sanitization = csi_sanitization.*ToMult;
    csi_Sanitization(:,:,i_t) = csi_sanitization;
end
relChannel_noSlope = reshape(csi_Sanitization, N_c, N_r, T);
sample_csi_trace_sanitized = relChannel_noSlope(:);

% MUSIC algorithm for estimating angle of arrival
% aoaEstimateMatrix is (nComps x 5) matrix where nComps is the number of paths in the environment. First column is ToF in ns and second column is AoA in degrees as defined in SpotFi paper
aoaEstimateMatrix = backscatterEstimationMusic(sample_csi_trace_sanitized, N_r, N_c, c, fc,...
                    T, f_gap, subcarrier_idx, d, paramRange, maxRapIters, useNoise, do_second_iter, ones(L+1));   
Dishat_spotFi(:,i_csi) = aoaEstimateMatrix(:,1); % ToF in nanoseconds
Thetahat_spotFi(:,i_csi) = aoaEstimateMatrix(:,2); % AoA in degrees
% Thetaerr_spotFi = abs(Thetahat_spotFi(end:-1:1).' - Theta)





end
Dishat_SpotFi = Dishat_spotFi(1:2,:);
Thetahat_SpotFi = Thetahat_spotFi(1:2,:);
figure;
plot(Dishat_SpotFi,Thetahat_SpotFi,'o')
axis([-40 40 -90 90])
xlabel('TOF','FontName','Times New Roman','FontSize',14);
ylabel('AOA','FontName','Times New Roman','FontSize',14);
set(gca,'xtick',-40:20:40,'ytick',-90:30:90)
grid on

DisTheta_Data = [Dishat_SpotFi(:) Thetahat_SpotFi(:)];
%% clustering

opts = statset('Display','iter');
[idx,C,sumd,d,midx,info] = kmedoids(DisTheta_Data,6,'Distance','mahalanobis','Options',opts);

figure;
plot(DisTheta_Data(idx==1,1),DisTheta_Data(idx==1,2),'r*','MarkerSize',10)
hold on
plot(DisTheta_Data(idx==2,1),DisTheta_Data(idx==2,2),'g*','MarkerSize',10)
hold on
plot(DisTheta_Data(idx==3,1),DisTheta_Data(idx==3,2),'b*','MarkerSize',10)
hold on
plot(DisTheta_Data(idx==4,1),DisTheta_Data(idx==4,2),'y*','MarkerSize',10)
hold on
plot(DisTheta_Data(idx==5,1),DisTheta_Data(idx==5,2),'c*','MarkerSize',10)
hold on
plot(DisTheta_Data(idx==6,1),DisTheta_Data(idx==6,2),'m*','MarkerSize',10)
hold on
plot(DisTheta_Data(idx==7,1),DisTheta_Data(idx==7,2),'k*','MarkerSize',10)
hold on
plot(C(:,1),C(:,2),'ko', 'MarkerSize',7,'LineWidth',1.5)
hold off

axis([-40 40 -90 90])
xlabel('TOF','FontName','Times New Roman','FontSize',14);
ylabel('AOA','FontName','Times New Roman','FontSize',14);
set(gca,'xtick',-40:20:40,'ytick',-90:30:90)
grid on


