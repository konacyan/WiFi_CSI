function TOF = MUSIC_tau(U, M, spacing, lambda )

c = 3e8; % speed of light
U = U./Column_wise_norm(U);
Proj = @(A)A*pinv(A'*A)*A';
% AOA estimation
Un = eye(M) - Proj(U);
% search


dis_spacing = 0.1;
dis_range = [1,15];
dis = dis_range(1):dis_spacing:dis_range(2); %distance search range, research step is 0.01m
Pmusic = zeros(1,length(dis)); % spectral


% !Steering vector function need to be modified here cause the antenna array
% !after smoothing is NOT equidistant among arrays



for i_dis = 1: length(dis)
	P_tau = exp(-1j*2*pi*[0:M-1]'*spacing *dis(i_dis)/c);
	Pmusic(i_dis) = 1./abs(P_tau'*Un*P_tau); % spectral

end

Pmusic=abs(Pmusic)./max(abs(Pmusic));

% % figure
% mesh(Pmusic);

% dis_idx = find(sum(Pmusic,1) == max(sum(Pmusic,1)));
% Pmusic_theta = Pmusic(:,dis_idx);
% figure;plot(theta,Pmusic_theta);


Pmusic = 10*log10(Pmusic);
[~, dis_idx] = find(Pmusic == max(Pmusic));

TOF = dis(dis_idx);




end



