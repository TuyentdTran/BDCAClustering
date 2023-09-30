% Find the matrix sum_{i=1}^m e_{r(i)*(x^{(r(i)}-w^{r(i)}:=Temp

function W = h1(X,C,R)
% C is an 2-by-m matrix containing centers of the target balls
% R is an 1-by-m array containing corresponding radius

[k,n]=size(X);
m=size(C,1);
W=zeros(k,n);

XA2 = zeros(m,k); % the j column, i row holds d(x_j,Omega_i)^2 %
ProjX = zeros(m,n,k); 

for i = 1:k
    ProjX(:,:,i) = projball2(X(i,:),C,R);
    XA2(:,i) =  vecnorm(ProjX(:,:,i) - X(i,:),2,2);
end
[~,I]=min(XA2,[],2);

for i=1:k   
    W(i,:) = sum(X(i,:) - ProjX(I==i,:,i),1);
end
end
  