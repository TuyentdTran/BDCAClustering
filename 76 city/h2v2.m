function W = h2v2( X , A )

[k,d]=size(X); m=size(A,1); 

W=zeros(k,d);
XA2norm = zeros(m,k); % the j column, i row holds \|x_j - a_i\|_2^2 %
for j = 1:d
    XA2norm = XA2norm + (A(:,j) - X(:,j)').^2;
end
[~,I] = min(XA2norm,[],2);
for i = 1:k % subtract the\sum_{i=1}^m e_{r(i)} \oprod (x^{r(i)} - a^i) term
    W(i,:) =  sum(X(i,:)- A(I==i,:),1); 
end
W = m*X - sum(A,1)-W; 
end
