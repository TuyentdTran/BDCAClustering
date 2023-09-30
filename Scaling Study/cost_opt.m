function c = cost_opt(A,X)

% An optimized cost penalty function that doesn't use DCA splitting to 
% evaluate 

%%
m=size(A,1);
k=size(X,1);
n=size(A,2);

XA2 = zeros(m,k); % the j column, i row holds \|x_j - a_i\|_2^2 %
for d = 1:n
    XA2 = XA2 + (A(:,d) - X(:,d)').^2;
end
c = 1/2*sum(min(XA2,[],2));
end