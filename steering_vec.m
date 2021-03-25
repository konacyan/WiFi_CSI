function a = steering_vec(theta,M,d,lambda)
% steering vector
% theta row vector
[m,n] = size(theta);
if m~=1
    disp('Error!theta needs to be a row vector')
    return  
end

a = exp(-1j*2*pi*[0:M-1]'*d*sin(theta)/lambda);
% a = a./norm(a);
end