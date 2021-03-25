function parameter = Improved_EF(data, paramRange, signal_space)  %%%data大小是data_rowxdata_col的矩阵

[data_row,data_col] = size(data);

%%%%%%%%40M带宽下的子载波索引%%%%%%%%%
% k = zeros(30,1);
% k(1:30,1) = -58: 4 : 58;
%%%%%%%%40M带宽下的子载波频率%%%%%%%%%
% f = zeros(30,1);
% f(15,1) = 5754.375;   %MHz
% f(16,1) = 5755.625;
% f(1:14,1) = 5736.875: 1.25 : 5753.125;
% f(17:30,1) = 5756.875: 1.25 : 5773.125;

%%%%%计算链路的平均多径数量%%%%%%%%%%%%%%
% for i = 1 : data_row
%     for j = 1 : data_col/30
%         path(i,j) = path_estimate(data(i,(j-1)*30+1:j*30));
% %         path(i,2) = path_estimate(data(i,31:60));
% %         path(i,3) = path_estimate(data(i,61:90));
%     end
% end
% path = ceil(mean(path)); %%%1x(data_col/30)行向量
%%%%%%%%%%%%%%%%%

% parameter = cell(1,data_col);
parameter = [];
%%%%%%%%对每一根天线分别进行滤波%%%%%%%%%%%%%%%
for i = 31 : 60
    parameter = [parameter ; EM1(data(:,i),0.0001,i-30,paramRange,2,signal_space)];
end
for i = 61 : 90
    parameter = [parameter ; EM1(data(:,i),0.0001,i-60,paramRange,3,signal_space)];
end
% for sub = 1 : data_col    %col个子载波要分别求解
% %     parameter{1,sub} = EM1(data(:,sub),0.01,path(1,ceil(sub/30)));
%     if mod(sub,30) == 0
% %        parameter{1,sub} = EM_new(data(:,sub),0.001,path(1,ceil(sub/30)),f(30,1),1.25,k(30,1));
%        parameter{1,sub} = EM1(data(:,sub),0.0001,path(1,ceil(sub/30)));
%        %parameter{1,sub} = EM2(data(:,sub),0.01,path(1,ceil(sub/30)));
%     else
% %        parameter{1,sub} = EM_new(data(:,sub),0.001,path(1,ceil(sub/30)),f(mod(sub,30),1),1.25,k(mod(sub,30),1));
%        parameter{1,sub} = EM1(data(:,sub),0.0001,path(1,ceil(sub/30)));
%        %parameter{1,sub} = EM2(data(:,sub),0.01,path(1,ceil(sub/30)));
%    end
% end
%%%%%%%%%%%%求得传播时间%%%%%%%%%%%%%%%
% parameter = ConverToTime(parameter);
%%%%%%%%%%%%挑选出LoS路径%%%%%%%%%%%%%%
% para = LoS(parameter);








%%%%%%从返回结果中找出直接路径%%%%%%%%%
%%%%%由于每条子载波的频率是一定的，而phi-e*pi*f*tau，所以最小值是直接路径%%%
%%%%%由于目前的相位值位于[-pi,pi]之间，对于小于0的，我们要加pi，使得其位于[0,2*pi]内，取最小
%%%%%同时，幅度值应大于0%%%%%%%%%%%%%%
% 
% for i = 1 : data_col
%     temp = parameter{1,i}(parameter{1,i}(:,1)>=0,:);
%     if isempty(temp)
%         amplitude(1,i) = 0;
%         temp1 = parameter{1,i};
%         [p_row,p_col] = size(parameter{1,i});
%         for j = 1 : p_row
%             if parameter{1,i}(j,2) < 0
%                 temp1(j,2) = temp1(j,2)+pi;
%             end
%         end
%         [minimum,index] = min(temp1(:,2));
%         phase(1,i) = parameter{1,i}(index,2);
%     else
%         temp1 = temp;
%         [t_row,t_col] = size(temp);
%         %%%%%将相位范围重置[0,2*pi]内%%%%%%%%%%
%         for j = 1 : t_row
%             if temp(j,2) < 0
%                 temp1(j,2) = temp1(j,2)+pi;
%             end
%         end
%         [minimum,index] = min(temp1(:,2));
%         amplitude(1,i) = temp(index,1);
%         phase(1,i) = temp(index,2);
%     end
% end