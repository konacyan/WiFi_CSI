function [S_T, N_T, Cortemplate, lags] = CorM(Tau, fs, K, I_1, beta, span, sps)
% y = rcosdesign(beta,span,sps);
% beta: roll-off factor [0,1]
% span: number of symbols [3]
% sps: samples per symbol, pulse duration [2.4ns]
% sps = ceil((fs/fb)/span);
% observation time [ 50 ]ns

y = raised_cosine_pulse(K, beta, span, sps, fs, 0);% signal template
ind = find(abs(y)>0);
S_T = zeros(K,K);
S_T(1,:) = y;
for k = 2:K-length(ind)+1
    % Positive k shifts toward the end of the dimension and negative K shifts toward the beginning.
    % dim = 1 to exchange rows, dim = 2 to exchange columns, and so on
    S_T(k,:) = circshift(y,k-1);
end
for k = K-length(ind)+2:K
    S_T(k,:) = [zeros(1,k-1) y(1:K-k+1)];
end
autocor = (y*S_T')'; 
ind = find(abs(autocor)>=1e-6);
Cortemplate = autocor(1:ind(end)); % correlation template(column)
lags = length(Cortemplate);


temp = fix(ceil(Tau(1)*fs)/I_1) - fix((ceil(Tau(1)*fs) + lags)/I_1);
if temp ~= 0
    Error('You are using an inappropriate I_1');
end

N_T = ceil(K/I_1)*I_1;
S_T = [S_T; zeros(N_T-K,K)];
% N_T = fix((K-length(ind)+1)/I_1)*I_1;


% figure
% hold on
% subplot(2,1,1)
% plot(1:60,y(1:60),'b-','linewidth',2);
% legend('Signal waveform (raised cosine)')
% set(gca,'FontName','Times New Roman', 'FontSize', 14,'xtick',[10:10:60],'ytick',[-0.25,0,0.25,0.5],'yticklabels',{'-0.2','0','0.25','0.5'}) 
% axis([1,60,-0.25,0.6])
% subplot(2,1,2)
% plot(1:100,autocor(1:100),'r-','linewidth',2);
% legend('Autocorrelation')
% set(gca,'FontName','Times New Roman', 'FontSize', 14,'xtick',[10:20:100])
% axis([1,100,-0.5,1.1])


end