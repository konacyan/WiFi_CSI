function y = tensorization(Y,I,mode)
% Y is matrix
% I is the size of the tensor
[M,N] = size(Y);
if cumprod(I) ~= M*N
    error('Unable to perform assignment because the size of the matrix is incompatible to the size of the tensor')
end

y = zeros(I);
if mode == 1
    for i=1:N
        y(:,:,i) = reshape(Y(:,i),I(1),I(2)); % column-wise
        % A = 1:10; B = reshape(A,[5,2]) [1:5 6:10]
    end
else
    for i=1:M
        y(:,:,i) = reshape(Y(i,:),I(1),I(2));
    end
end
    

end