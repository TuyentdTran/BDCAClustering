function [c,I]=cost(A,X,projfun)
%%
%%
m=size(A,1);
k=size(X,1);
% n=size(A,2);

XA2 = zeros(m,k); % the j column, i row holds d(x_j - Omega_i)^2 %
for i = 1:k
    XA2(:,i) =  vecnorm(projfun(X(i,:)) - X(i,:),2,2).^2;
end

[~,I]=min(XA2,[],2);
c=sum(min(XA2,[],2));
end