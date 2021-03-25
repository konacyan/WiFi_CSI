%%%%%%%%%%����һԪ���Իع飬������������SFO������%%%%%%%%%%%%
%%%%%%%%%%���ߣ��Բ���%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%�Ժ�����ʽ����ʱdata��1x30��������%%%%%%%%%%%%%%%%%

function phase = phase_calibration(data) %�����������λ����ֵ

% data = xlsread('/Users/ucaszbb/Desktop/20191227phase.xlsx','Sheet1','A1:AD3');
[row,col] = size(data);

%%%%%�����%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1 : row
    data(i,1:30) = unwrap(data(i,1:30));
end

%%%%%��������ƺ����λ%%%%%%%%%%%%%%%%
% plot(1:30,data(1,1:30),'-ro','linewidth',2,'markeredgecolor','b','markerfacecolor','m','markersize',10);
% hold on;
% plot(1:30,data(2,1:30),'-go','linewidth',2,'markeredgecolor','b','markerfacecolor','y','markersize',10);
% hold on;
% plot(1:30,data(3,1:30),'-bo','linewidth',2,'markeredgecolor','b','markerfacecolor','k','markersize',10);


%%%%%%%%20M�����µ����ز�����%%%%%%%%%
% k = zeros(30,1);
% k(1:14,1) = -28: 2 : -2;
% k(15,1) = -1;
% k(16:29,1) = 1 : 2 : 27;
% k(30,1) = 28;
%%%%%%%%40M�����µ����ز�����%%%%%%%%%
k = zeros(30,1);
k(1:30,1) = -58: 4 : 58;

phase = zeros(row,col);    %��֤����30
%phase = zeros(3,30);

for i = 1 : row   
    a = (data(i,30)-data(i,1))/(k(30,1)-k(1,1));
    b = sum(data(i,1:30),2)/30;
    for j = 1 : 30   
        phase(i,j) = data(i,j)-(a*k(j,1)+b);
    end
end
% figure();
% for i = 1 : row
%     plot(1:30,phase(i,1:30),'-o','linewidth',2,'markersize',10);
%     hold on;
% end

% plot(1:30,phase(1,1:30),'-ro','linewidth',2,'markeredgecolor','b','markerfacecolor','m','markersize',10);
% hold on;
% plot(1:30,phase(2,1:30),'-go','linewidth',2,'markeredgecolor','b','markerfacecolor','y','markersize',10);
% hold on;
% plot(1:30,phase(3,1:30),'-bo','linewidth',2,'markeredgecolor','b','markerfacecolor','k','markersize',10);

% % legend('RX=1,packet=1ԭʼ���������','RX=1,packet=2ԭʼ���������','RX=1,packet=3ԭʼ���������','RX=1,packet=1ƫ������������','RX=1,packet=2ƫ������������','RX=1,packet=3ƫ������������');
% legend({'RX=1,packet=1','RX=1,packet=2','RX=1,packet=3'},'FontSize',15);
% xlabel('subcarrier index','FontSize',18);
% ylabel('unwrapped phase (radian)','FontSize',18);
% axis([0,30,-4,4]);