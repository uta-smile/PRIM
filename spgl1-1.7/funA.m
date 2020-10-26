function out = funA(X,Phi,PF,T,m,n)
N=m*n;
for t=1:T
    x1 = X(:,(t-1)*2+1);
    x2 = X(:,t*2);
    x = x1+1i*x2;
    x = reshape(x,m,n);
    y = PF{t}*(Phi'*x);
    y = y(:);
    out(:,(t-1)*2+1) = real(y);
    out(:,t*2) = imag(y);
end