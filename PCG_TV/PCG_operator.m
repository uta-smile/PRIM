function [c,k,res,err] = PCG_operator(A,b,invP,kmax,c,tol,prec,x0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The revised PCG solver. A can be both an operator or a numerical matrix.
%%% Solve Ax = b with preconditioner P
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



nb = norm(b,'fro'); 
res = zeros(kmax+1,1);
err = res;

N = length(b);

r = b - A*c;                       % r = b-A*c;
res(1) = norm(r,'fro')^2;

if prec==0                                        % proceed without preconditioner
  z  = r;
  for k=1:kmax
    w = A*z;
    alpha = res(k)/(z(:)'*w(:));
    c = c + alpha*z;
    r = r - alpha*w;
    res(k+1) = norm(r,'fro')^2;
    beta = res(k+1)/res(k);
    z = r + beta*z;
    if sqrt(res(k+1))<tol*nb
      break
    end
    err(k)=norm(c-x0)/norm(x0);
%       err(k)=norm(c-x0)/length(x0);
  end
end  %0
  
 if prec==1                                                      % use preconditioner
  z = invP*r;
  rs = r(:)'*z(:);

  for k=1:kmax
    w = A*z;
    alpha = rs/(z(:)'*w(:));
    c = c + alpha*z;
    r = r - alpha*w;
    res(k+1) = norm(r,'fro')^2;
    s = invP*r;
   % s=S\r;
   % ds =norm(s -s1);
    rs_old = rs;
    rs = r(:)'*s(:);
    beta = rs/rs_old;
    z = s + beta*z;
    if sqrt(res(k+1))<tol*nb
      break
    end
    err(k)=norm(c-x0)/norm(x0);
%       err(k)=norm(c-x0)/length(x0);
  end
  
 end%2
    
end   % function
