function AOA = MUSIC_joint(U, n_c,spacing, M, m, d, lambda )
n_c = n_c +1;% virtaul array(n_c) add physical array(1)
c = 3e8; % speed of light
U = U./Column_wise_norm(U);
Proj = @(A)A*pinv(A'*A)*A';
% AOA estimation
Un = eye(M) - Proj(U);
% search


theta_spacing = 0.5;
theta_range = [-90,90];
dis_spacing = 2;
dis_range = [500,1500];
theta = (theta_range(1) + dis_spacing):theta_spacing:(theta_range(2) - dis_spacing); % angle search range, research step is 0.01 degree
dis = dis_range(1):dis_spacing:dis_range(2); %distance search range, research step is 0.01m
Pmusic = zeros(length(theta),length(dis)); % spectral


% !Steering vector function need to be modified here cause the antenna array
% !after smoothing is NOT equidistant among arrays

for i_theta=1:length(theta)
    phy_array = steering_vec(theta(i_theta)*pi/180,m,d,lambda); %physical array
%     phy_array = exp(-1j*2*pi*[0:m-1]'*d*sin(i_theta)/lambda); %physical array
    %         Pmusic{l}(i) = 1./abs(temp'*(Un*Un')*temp); % spectral
    for i_dis = 1: length(dis)
        vir_array = exp(-1j*2*pi*[0:n_c-1]'*spacing *dis(i_dis)/c);
        for i_thelen = 1: length(phy_array)
            tol_array((i_thelen-1)*length(vir_array)+1:i_thelen*length(vir_array),1) = phy_array(i_thelen).*vir_array;
        end
        Pmusic(i_theta,i_dis) = 1./abs(tol_array'*Un*tol_array); % spectral
    end
    
    
end

Pmusic=abs(Pmusic)./max(max(abs(Pmusic)));
%
% figure
plot(Pmusic)
Pmusic = 10*log10(Pmusic);
[theta_idx, dis_idx] = find(Pmusic == max(max(Pmusic)));
% theta_est = theta_range(1) + theta_idx * theta_spacing;
AOA = theta(theta_idx);
TOF = dis(dis_idx);
% dis_idx = find(sum(Pmusic,1) == max(sum(Pmusic,1)));

end



