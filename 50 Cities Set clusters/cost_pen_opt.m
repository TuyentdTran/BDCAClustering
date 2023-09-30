function c=cost_pen_opt(C,R,X,tau,proj_c)

% An optimized cost penalty function that doesn't use DCA splitting to 
% evaluate 

%%
m=size(C,1);
k=size(X,1);


XA2 = zeros(m,k); % the j column, i row holds d(x_j - Omega_i)^2 %
for i = 1:k
    XA2(:,i) =  vecnorm(projball2(X(i,:),C,R) - X(i,:),2,2).^2;
end
%for general 2 constraints consecutive

temp = proj_c(X);
ch_2 = 0;
for i=1:k
    ch_2 = ch_2 +  sum(vecnorm(X(i,:)-temp(2*i-1:2*i,:),2,2).^2); 
end

c = 1/2*sum(min(XA2,[],2))+tau/2*ch_2;

end