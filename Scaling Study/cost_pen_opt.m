function c=cost_pen_opt(A,X,tau,q,projfun)

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

%hard coded for the given problem. Namely this assumes a constraint has a
%single domain associated with it.
c = 1/2*sum(min(XA2,[],2))+tau/2*(sum(vecnorm(X-projfun(X),2,2).^2)); % this 2nd term is hardcoded
end