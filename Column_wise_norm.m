function y = Column_wise_norm(X)

[~,N] = size(X);
y = zeros(1,N);
for i=1:N
   y(i) = norm(X(:,i)); 
end
end