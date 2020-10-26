function out = funPhiT(X,Phi,T,m,n)
N=m*n;
X = reshape(X,N,T*2);
for t=1:T
    x1 = X(:,(t-1)*2+1);
    x2 = X(:,t*2);
    x = x1+1i*x2;
    x = reshape(x,m,n);
    y = Phi'*x;
    out(:,:,t) = y;
    
end