clc
clear 
addpath('tensorlab_2016-03-28');
addpath('linux-80211n-csitool-supplementary');
%% code parameter
simulation_flag = 0;
plot_flag = 0;
phaseCalibration_flag = 0;
%% system setting
% fc = 2.472e9; % Channel 13
fc = 5.32e9; % Channel 64
f_gap = 312.5e3;% subcarrier spacing
c = 3e8; % speed of light
lambda = c/fc;  % wavelength(m)
d = 5.35e-2./2;
A = @(theta,M)steering_vec(theta/180*pi,M,d,lambda);

% subcarrier_idx = [-28:2:-2,-1,1:2:27,28];% 20MHT mode, carriers index in IEEE 802.11n
subcarrier_idx = -58:4:58;% 40MHT mode, carriers index in IEEE 802.11n
subcarrier_spacing = 4 * f_gap;
f_carrier = fc + f_gap.*subcarrier_idx;

%% load data
% filepath = 'Data/0310_B/A2B0003_0.dat';
filepath = 'Data/0301_B/log0301A2B11_0.dat';
% filepath = 'Data/0316_AP/0316log23_0.dat';
WLAN = wifi_config();
[csi_data,csi_trace] = gene_csi(filepath);
% save('log0129A2B64HT20-11_0.mat', 'csi_trace');
% save('test1_90','csi_data')

N_pack = 20; % package number
for i_csi = 1:N_pack
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
    csi_SNR_rxC = zeros(N_r,N_c);
    for j=1:N_r   
        for k=1:N_c    
            csi_SNR_rxA(j,k)=csi(1,j,k);
            csi_SNR_rxB(j,k)=csi(2,j,k);
            csi_SNR_rxC(j,k)=csi(3,j,k);
        end
    end
    figure;
    subplot(1,3,1)
    plot(db(abs(csi_SNR_rxA.')))
    legend('RX Antenna A', 'RX Antenna B', 'RX Antenna C');
    xlabel('Subcarrier index');
    ylabel('SNR [dB]')
    subplot(1,3,2)
    plot(db(abs(csi_SNR_rxB.')))
    legend('RX Antenna A', 'RX Antenna B', 'RX Antenna C');
    xlabel('Subcarrier index');
    ylabel('SNR [dB]')
    subplot(1,3,3)
    plot(db(abs(csi_SNR_rxC.')))
    legend('RX Antenna A', 'RX Antenna B', 'RX Antenna C');
    xlabel('Subcarrier index');
    ylabel('SNR [dB]')

% plot unwrapped phase 
    N_pack = 1;% first package
    for i_t = 1:N_t
        for i_r = 1:N_r
            csi_pack = csi_data{N_pack};
            csi_subcarriers(:,((i_t - 1)*N_t + i_r)) = squeeze(csi_pack(i_t,i_r,:)).';
        end
    end
    csi_phase=angle(csi_subcarriers);
    csi_unwrapphase = unwrap(csi_phase);
    figure;
    plot(180*csi_phase/pi,'.','Markersize',25);
    xlabel('Subcarrier Index');
    ylabel('Phase [deg]')

% plot unwrap phase
    figure;
    for i_s = 1:N_r*N_t
            subplot(N_t, N_r,i_s)
            plot(180/pi*csi_unwrapphase(:,i_s));
            xlabel('Subcarrier Index');
            ylabel('Unwrap Phase');
            tx_num = fix((i_s-1)/N_t)+1;
            rx_num = mod((i_s-1),N_t)+1;
            title(['TX:',num2str(tx_num) '  RX:',num2str(rx_num)])
    end
 
end

%% phase calibration
if phaseCalibration_flag == 1
    figure;
    coefficient = zeros(N_t*N_r,2);
    csi_poly = zeros(N_t*N_r,N_c);
    for i_t = 1: N_t*N_r

        coefficient(i_t,:) = polyfit(1:N_c,csi_unwrapphase(:,i_t).',1);
        csi_poly(i_t,:) = polyval(coefficient(i_t,:),1:N_c);
        subplot(N_t,N_r,i_t)
        plot(1:N_c,csi_unwrapphase(:,i_t),'o',1:N_c,csi_poly(i_t,:),'-');

    end

    csi_modified_phase = csi_unwrapphase - csi_poly.';
    csi_STO = -1*coefficient(:,1)./(2*pi*subcarrier_spacing);% can be further used to get the multipath number(SpotFi, page 8)
    figure;
    for i_t = 1: N_t*N_r
        subplot(N_t,N_r,i_t)
        plot(1:N_c,csi_modified_phase(:,i_t));
    end


    csi_amp = abs(csi_subcarriers);
    csi_modified = csi_amp.*(cos(csi_modified_phase) + 1i*sin(csi_modified_phase));
    csi_Modified = zeros(N_t, N_r, N_c);
    for i_t = 1: N_t
        for i_r = 1:N_r
            for i_c = 1:N_c
                csi_Modified(i_t,i_r,i_c) = csi_modified(i_c,(i_t-1)*N_t+i_r);
            end
        end
    end
    csi_unModified = csi;
    csi = csi_Modified;
end

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
L = 2;

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
                    T, f_gap, subcarrier_idx, d, paramRange, maxRapIters, useNoise, do_second_iter, ones(L+1)) ;  
if size(aoaEstimateMatrix,1) ~= L
    aoaestimateMatrix = zeros(L,5);
    matrixlength = size(aoaEstimateMatrix,1);
    for i_a = 1:matrixlength
        aoaestimateMatrix(i_a,:) = aoaEstimateMatrix(i_a,:);
    end
    aoaEstimateMatrix = aoaestimateMatrix;
end
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


