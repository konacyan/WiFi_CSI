function WLAN=wifi_config()
WLAN.freq = 5.745 * 10^9; %信道频率
WLAN.sub_freq_delta = 312500*4;%子载波频率间隔
WLAN.array=zeros(3,3);
% R=WLAN.ant_dis*sqrt(3)/2;
% WLAN.array(1,:)=[R 0 0];
% WLAN.array(2,:)=[-R/2 sqrt(3)*R/2 0];
% WLAN.array(3,:)=[-R/2 -sqrt(3)*R/2 0];
WLAN.array(1,:)=[-0.025 0 0];
WLAN.array(2,:)=[0 0 0];
WLAN.array(3,:)=[0.025 0 0];
end