function parameter = Improved_EF(data, paramRange, signal_space)  %%%data��С��data_rowxdata_col�ľ���

[data_row,data_col] = size(data);

%%%%%%%%40M�����µ����ز�����%%%%%%%%%
% k = zeros(30,1);
% k(1:30,1) = -58: 4 : 58;
%%%%%%%%40M�����µ����ز�Ƶ��%%%%%%%%%
% f = zeros(30,1);
% f(15,1) = 5754.375;   %MHz
% f(16,1) = 5755.625;
% f(1:14,1) = 5736.875: 1.25 : 5753.125;
% f(17:30,1) = 5756.875: 1.25 : 5773.125;

%%%%%������·��ƽ���ྶ����%%%%%%%%%%%%%%
% for i = 1 : data_row
%     for j = 1 : data_col/30
%         path(i,j) = path_estimate(data(i,(j-1)*30+1:j*30));
% %         path(i,2) = path_estimate(data(i,31:60));
% %         path(i,3) = path_estimate(data(i,61:90));
%     end
% end
% path = ceil(mean(path)); %%%1x(data_col/30)������
%%%%%%%%%%%%%%%%%

% parameter = cell(1,data_col);
parameter = [];
%%%%%%%%��ÿһ�����߷ֱ�����˲�%%%%%%%%%%%%%%%
for i = 31 : 60
    parameter = [parameter ; EM1(data(:,i),0.0001,i-30,paramRange,2,signal_space)];
end
for i = 61 : 90
    parameter = [parameter ; EM1(data(:,i),0.0001,i-60,paramRange,3,signal_space)];
end
% for sub = 1 : data_col    %col�����ز�Ҫ�ֱ����
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
%%%%%%%%%%%%��ô���ʱ��%%%%%%%%%%%%%%%
% parameter = ConverToTime(parameter);
%%%%%%%%%%%%��ѡ��LoS·��%%%%%%%%%%%%%%
% para = LoS(parameter);








%%%%%%�ӷ��ؽ�����ҳ�ֱ��·��%%%%%%%%%
%%%%%����ÿ�����ز���Ƶ����һ���ģ���phi-e*pi*f*tau��������Сֵ��ֱ��·��%%%
%%%%%����Ŀǰ����λֵλ��[-pi,pi]֮�䣬����С��0�ģ�����Ҫ��pi��ʹ����λ��[0,2*pi]�ڣ�ȡ��С
%%%%%ͬʱ������ֵӦ����0%%%%%%%%%%%%%%
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
%         %%%%%����λ��Χ����[0,2*pi]��%%%%%%%%%%
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