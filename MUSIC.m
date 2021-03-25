function AOA = MUSIC(U, M, d, lambda )


U = U./Column_wise_norm(U);
Proj = @(A)A*pinv(A'*A)*A';
% AOA estimation
Un = eye(M) - Proj(U);
% search
% the = theta-Dtheta:0.01:theta + Dtheta; % angle search range
the = -89.99:0.01:89.99; % angle search range, research step is 0.01 degree
Pmusic = zeros(length(the),1); % spectral


% !Steering vector function need to be modified here cause the antenna array
% !after smoothing is NOT equidistant among arrays

for i=1:length(the)
    temp = steering_vec(the(i)*pi/180,M,d,lambda);
    %         Pmusic{l}(i) = 1./abs(temp'*(Un*Un')*temp); % spectral
    Pmusic(i) = 1./abs(temp'*Un*temp); % spectral
end

Pmusic=abs(Pmusic);
%
% figure
plot(the,Pmusic)
Pmusic = 10*log10(Pmusic/max(Pmusic));
[~, ind] = max(Pmusic);
AOA = the(ind);


end



