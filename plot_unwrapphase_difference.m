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
Lambda = c./f_carrier;

%% load data

% filenum = [11 12 13 20 21 22 23 -10-11 -12 -13 -20 -21 -22 -23 -30 -31 -32 -33];
% for i_file =1: length(filenum)
%     fileNum = num2str(filenum(i_file));

filename = '0316log-10_1';
filepath = ['Data/0316_AP/',filename,'.dat'];
WLAN = wifi_config();
[csi_data,csi_trace] = gene_csi(filepath);
% save([filename,'.mat'], 'csi_trace');

% Groundtruth
groundtruth_x = str2num(filename(end-3));
groundtruth_y = str2num(filename(end-2));
angle_true =  90 - 180/pi*atan(groundtruth_y/groundtruth_x)
phase_difference_true = 2*pi*d/lambda*sin(angle_true*pi/180)
% only the first data is used

% for i_csi = 1:length(csi_data)
%     csi_size =  size(csi_data{i_csi});
%     if csi_size(1,1) > 1
%         csi = csi_data{i_csi};
%         csi_detail = csi_trace{i_csi};
%         break;
%     end
% end
csi_packnum = length(csi_data);
for i_csi = 1:20

csi = csi_data{i_csi};
csi_detail = csi_trace{i_csi};

% sort csi
[N_t, N_r, N_c] = size(csi);
csi_perm = csi_detail.perm;
csi_perm_idx = zeros(1,N_r);
for i_r = 1:N_r
    csi_perm_idx(i_r) = find(csi_detail.perm == i_r);
end
csi_sort = zeros(N_t,N_r,N_c);

for i_r = 1:N_r
    csi_sort(:,i_r,:) = csi(:,csi_perm_idx(i_r),:);
end
csi = csi_sort;
% if csi_detail.perm == [1 3 2]
%     csi_sort = zeros(N_t,N_r,N_c);
%     csi_sort(:,1,:)=csi(:,1,:);
%     csi_sort(:,3,:)=csi(:,2,:);
%     csi_sort(:,2,:)=csi(:,3,:);
%     csi = csi_sort;
% end

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
%     N_pack = 1;% first package
    for i_t = 1:N_t
        for i_r = 1:N_r
            csi_pack = csi;
            csi_subcarriers(:,((i_t - 1)*N_r + i_r)) = squeeze(csi_pack(i_t,i_r,:)).';
        end
    end
    csi_amp = abs(csi_subcarriers);
    csi_normal = csi_subcarriers./csi_amp;
    csi_phase=angle(csi_normal);
    csi_unwrapphase = unwrap(csi_phase);
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
    

    
    
    
%     [~,csi_unwrapphase_idx] = sort(sum(csi_unwrapphase(:,1:3),1));
    csi_Difference = zeros(N_c,N_t*N_r);
    for i_t = 1:N_t
        csi_Difference(:,(i_t - 1)*N_r + 1)= csi_unwrapphase(:,(i_t - 1)*N_r + 1) - csi_unwrapphase(:,(i_t - 1)*N_r + 2);
        csi_Difference(:,(i_t - 1)*N_r + 2)= csi_unwrapphase(:,(i_t - 1)*N_r + 1) - csi_unwrapphase(:,(i_t - 1)*N_r + 3); 
        csi_Difference(:,(i_t - 1)*N_r + 3)= csi_unwrapphase(:,(i_t - 1)*N_r + 2) - csi_unwrapphase(:,(i_t - 1)*N_r + 3); 
    end
    
    
    if mean(csi_Difference(:,4)) < 0 &&  mean(csi_Difference(:,6)) > 0 
        csi_unwrapphase(:,5) = csi_unwrapphase(:,5) - 2*pi;
        csi_Difference(:,4)= csi_unwrapphase(:,4) - csi_unwrapphase(:,5);
        csi_Difference(:,6)= csi_unwrapphase(:,5) - csi_unwrapphase(:,6);
    
    elseif mean(csi_Difference(:,1)) < 0 &&  mean(csi_Difference(:,3)) > 0 
        csi_unwrapphase(:,2) = csi_unwrapphase(:,2) - 2*pi;
        csi_Difference(:,1)= csi_unwrapphase(:,1) - csi_unwrapphase(:,2);
        csi_Difference(:,3)= csi_unwrapphase(:,2) - csi_unwrapphase(:,3);
        
    end
    


    % calculate the true value
    phase_diff_true = 2*pi*d*sin(angle_true/180*pi)./Lambda;
    phase_diff2_true = 2*pi*2*d*sin(angle_true/180*pi)./Lambda;
    coefficient_true = polyfit(subcarrier_idx,phase_diff_true,1);
    coefficient2_true = polyfit(subcarrier_idx,phase_diff2_true,1);
    ft1 = fittype(@(b,x) coefficient_true(1)*x+b );
    ft2 = fittype(@(a,x) coefficient2_true(1)*x+a );
    for i_ft = 1:N_t
        fit_1_2 = fit(subcarrier_idx.',csi_Difference(:,(i_ft - 1)*N_r + 1),ft1,'start',0);
        coefficient(i_ft,1) = fit_1_2(0);
        fit_2_3 = fit(subcarrier_idx.',csi_Difference(:,(i_ft - 1)*N_r + 3),ft1,'start',0);
        coefficient(i_ft,2) = fit_2_3(0);
        fit_1_3 = fit(subcarrier_idx.',csi_Difference(:,(i_ft - 1)*N_r + 2),ft2,'start',0);
        coefficient2(i_ft,1) = fit_1_3(0);
    end

    % eliminate gap
    coefficient_gap = coefficient_true(2) - coefficient;
    coefficient2_gap = coefficient2_true(2) - coefficient2;
    Csi_unwrapphase = zeros(N_c, N_t*N_r);
    for i_t = 1:N_t
        Csi_unwrapphase(:,(i_t - 1)*N_r + 1)= csi_unwrapphase(:,(i_t - 1)*N_r + 1); 
        Csi_unwrapphase(:,(i_t - 1)*N_r + 2)= csi_unwrapphase(:,(i_t - 1)*N_r + 2) - coefficient_gap(i_t, 1); 
        Csi_unwrapphase(:,(i_t - 1)*N_r + 3)= csi_unwrapphase(:,(i_t - 1)*N_r + 3) - coefficient2_gap(i_t, 1); 
    end

    %calculate the new difference
    Csi_Difference = zeros(N_c,N_t*N_r);
    for i_t = 1:N_t
        Csi_Difference(:,(i_t - 1)*N_r + 1)= Csi_unwrapphase(:,(i_t - 1)*N_r + 1) - Csi_unwrapphase(:,(i_t - 1)*N_r + 2);
        Csi_Difference(:,(i_t - 1)*N_r + 2)= Csi_unwrapphase(:,(i_t - 1)*N_r + 1) - Csi_unwrapphase(:,(i_t - 1)*N_r + 3); 
        Csi_Difference(:,(i_t - 1)*N_r + 3)= Csi_unwrapphase(:,(i_t - 1)*N_r + 2) - Csi_unwrapphase(:,(i_t - 1)*N_r + 3); 
    end
        % plot the difference
    figure(1);
    set(figure(1),'Position',[100,100,800,400]);
    for i_t = 1:N_t
            subplot(1,2,i_t)
            plot(subcarrier_idx, Csi_unwrapphase(:,(i_t - 1)*N_r + 1),'r*-','linewidth',1.2);
            hold on;
            plot(subcarrier_idx, Csi_unwrapphase(:,(i_t - 1)*N_r + 2),'bo-','linewidth',1.2);
            hold on;
            plot(subcarrier_idx, Csi_unwrapphase(:,(i_t - 1)*N_r + 3),'k+-','linewidth',1.2);  
            hold off;
            xlabel('Subcarrier Index','FontName','Times New Roman','FontSize',14);
            ylabel('Unwrap Phase','FontName','Times New Roman','FontSize',14);
            l1 = legend('RX 1','RX 2','RX 3');
            set(l1, 'FontName','Times New Roman');
            axis([-28 28 -25 5]);
            title(['Pack Num: ',num2str(i_csi)],'FontName','Times New Roman','FontSize',14);
    end
    % plot gif
    [fig_A,fig_mapA] = rgb2ind(frame2im(getframe((gcf))),256);
    if i_csi == 1
        imwrite(fig_A,fig_mapA,'unwrapphase.gif','gif','LoopCount',Inf,'DelayTime',0.2);
    else 
        imwrite(fig_A,fig_mapA,'unwrapphase.gif','gif','WriteMode','append','DelayTime',0.2);
    end
    
    
    figure(2);
    set(figure(2),'Position',[50,50,800,400]);
    for i_t = 1:N_t
            subplot(1,2,i_t)
            plot(subcarrier_idx, Csi_Difference(:,(i_t - 1)*N_r + 1),'r*-','linewidth',1.2);
            hold on;
            plot(subcarrier_idx, Csi_Difference(:,(i_t - 1)*N_r + 2),'bo-','linewidth',1.2);
            hold on;
            plot(subcarrier_idx, Csi_Difference(:,(i_t - 1)*N_r + 3),'k+-','linewidth',1.2);  
            hold off;
            xlabel('Subcarrier Index','FontName','Times New Roman','FontSize',14);
            ylabel('Unwrap Phase','FontName','Times New Roman','FontSize',14);
            l2 = legend('RX 1 - RX 2','RX 1 - RX 3','RX 2 - RX 3');
            set(l2, 'FontName','Times New Roman');
            axis([-28 28 1 8]);
            title(['Pack Num: ',num2str(i_csi)],'FontName','Times New Roman','FontSize',14);
            
    end
    
    [fig_B,fig_mapB] = rgb2ind(frame2im(getframe((gcf))),256);
    if i_csi == 1
        imwrite(fig_B,fig_mapB,'phasediff.gif','gif','LoopCount',Inf,'DelayTime',0.2);
    else 
        imwrite(fig_B,fig_mapB,'phasediff.gif','gif','WriteMode','append','DelayTime',0.2);
    end
    
end

end