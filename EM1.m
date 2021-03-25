function parameter = EM1(data,threshold,sub,paramRange,antenna_num,path_number)
%%%%%%%%%%%%%%%%%%%
%parameter 存储返回的估计参数
%data data_colx1测量值，列向量
%threshold 阈值，小于阈值停止下一轮计算
%%%%%%%%%%%%%%%%%%%

%%%%%%对待估计参数有一个初始赋值%%%%%%%
% alpha = mean(abs(data),1)*ones(path_number,1);
alpha =rand(path_number,1)*10;
tau = rand(path_number,1)*10^-8;  %%%注意, 传播时间
theta = rand(path_number,1)*100;

alpha_ini = alpha;
tau_ini = tau;
theta_ini = theta;

x = zeros(path_number,length(data));
P = zeros(path_number,1);

while true
    for i = 1 : path_number
        %%%%%%%%E步骤%%%%%%%%
        P = alpha .* exp(-1i*2*pi*(paramRange.frequency_subcarrier(sub,1).*tau+(paramRange.frequency*(antenna_num-1)*paramRange.antenna_space.*cos(theta*pi/180))/paramRange.speed_light));
%         for k = 1 : path_number
%             P(k,1) = alpha(k,1) * (cos(phase(k,1))+1i*sin(phase(k,1)));
%         end
        for j = 1 : length(data) %列
            x(:,j) = P + (data(j,1)-sum(P)) * (1/path_number)*ones(path_number,1);      
        end
        %%%%%%%%%M步骤%%%%%%%%%%%%%%
        tau_ini(i,1) = tau(i,1);
        alpha_ini(i,1) = alpha(i,1);
        theta_ini(i,1) = theta(i,1);
        
        temp_part1 = imag(sum(conj(x(i,:)),2)/sum(x(i,:),2));
        temp_part2 = real(sum(conj(x(i,:)),2)/sum(x(i,:),2));
        C = atan(temp_part1/temp_part2)/(4*pi);
        
        
        tau(i,1) = (C-paramRange.frequency*((antenna_num-1)*paramRange.antenna_space*cos(theta(i,1)*pi/180)/paramRange.speed_light))/paramRange.frequency_subcarrier(sub,1);
        theta(i,1) = acosd(((C-paramRange.frequency_subcarrier(sub,1)*tau(i,1))*paramRange.speed_light)/(paramRange.frequency*(antenna_num-1)*paramRange.antenna_space));
        
        temp_part3 = sum(conj(x(i,:)),2);
        temp_part4 = sum(x(i,:),2);
        temp_part5 = exp(-1i*2*pi*((paramRange.frequency_subcarrier(sub,1)*tau(i,1)+(paramRange.frequency*(antenna_num-1)*paramRange.antenna_space*cos(theta(i,1)*pi/180))/paramRange.speed_light)));
        temp_part6 = exp(1i*2*pi*((paramRange.frequency_subcarrier(sub,1)*tau(i,1)+(paramRange.frequency*(antenna_num-1)*paramRange.antenna_space*cos(theta(i,1)*pi/180))/paramRange.speed_light)));
        
        alpha(i,1) = (0.5/length(data))*(temp_part3*temp_part5+temp_part4*temp_part6);
    end
    
    if pdist([tau_ini';tau'],'euclidean') < 10^-8 && pdist([alpha_ini';alpha'],'euclidean') < 0.1 && pdist([theta_ini';theta'],'euclidean') < 0.1
        break;
    end
end



parameter = [alpha, tau, theta];