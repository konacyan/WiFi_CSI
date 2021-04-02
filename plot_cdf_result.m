clc;
close all;
clear all;
%% Load Data
filename = 'MultipackResult03_2';
load(['Result/',filename,'.mat']);
Data_length = 400;
Theta_result_Tensor = Theta_result_tensor(1:400);
Theta_result_SpotFi = min(abs(Theta_result_spotFi(:,1:400)));



%% Plot CDF Figure
figure;
cdf_tensor = cdfplot(abs(Theta_result_Tensor));
set(cdf_tensor,'LineStyle','-','Color','b','Linewidth',1.5);
hold on;
cdf_spotfi = cdfplot(abs(Theta_result_SpotFi));
set(cdf_spotfi,'LineStyle','-','Color','r','Linewidth',1.5);
hold off;
xlabel('Error in degrees','FontName','Times New Roman','FontSize',16);
ylabel('Empirical   CDF','FontName','Times New Roman','FontSize',16);
l1 = legend('Tensor','SpotFi');
set(l1, 'FontName','Times New Roman','FontSize',14);
set(gca,'xtick',0:10:50,'xticklabel',{'0','10','20','30','40','50'},'FontName','Times New Roman','FontSize',14);
set(gca,'ytick',0:0.2:1,'yticklabel',{'0','0.2','0.4','0.6','0.8','1'},'FontName','Times New Roman','FontSize',14);
axis([0 50 0 1]);
grid on;
title(' ');
saveas(gcf,['Figure/',filename,'.jpg'])
saveas(gcf,['Figure/',filename,'.fig'])

