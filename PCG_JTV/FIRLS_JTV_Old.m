function output = FIRLS_JTV(input,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
l=input.l; u=input.u;
if((l==-Inf)&&(u==Inf))
    project=@(x)x;
elseif (isfinite(l)&&(u==Inf))
    project=@(x)(((l<x).*x)+(l*(x<=l)));
elseif (isfinite(u)&&(l==-Inf))
    project=@(x)(((x<u).*x)+((x>=u)*u));
elseif ((isfinite(u)&&isfinite(l))&&(l<u))
    project=@(x)(((l<x)&(x<u)).*x)+((x>=u)*u)+(l*(x<=l));
else
    error('lower and upper bound l,u should satisfy l<u');
end
ratio =  input.ratio;

n1=input.n1; n2=input.n2; N=n1*n2;
A=input.A; b=input.b;
T = size(b,2);
% At = input.At;
% if T==1
%     Atb = A'*b;
% else
    for t=1:T
        Atb(:,t)=A{t}'*b{t};
    end
% end
% l = n1;
f = 0;
y= Atb;
alpha = input.alpha;
iterno=0;
if isnumeric(A{t})
    for t=1:T
        AtA{t}=A{t}'*A{t};
    end
else
    for t=1:T
        AtA{t} =  A_operator(@(x) A{t}'*(A{t}*x), @(x) A{t}'*(A{t}*x));
    end
end
Q1 = input.Q1;
Q2 = input.Q2;
In = speye(N);

t00 = cputime; 
eps = 1e-14;
for itr = 1:input.no    % total iter counter
    iterno=iterno+1;
    
    yp=y;
    D1 = zeros(size(y(:,1)));
    D2 = D1;
    for t =1:T
        D1=D1 + (abs(Q1*y(:,t))+eps).^(2);
        D2=D2 + (abs(Q2*y(:,t))+eps).^(2);
    end
    D = (D1+D2).^(-0.5);
    
    P=ratio*In;
    tp=[(1:N)' (1:N)' D];
    dd=spconvert(tp);
    tp =  alpha*((Q1'*dd*Q1)+(Q2'*dd*Q2));
    
    
    for t=1:T
        S = A_operator(@(x) AtA{t}*x + tp*x, @(x) AtA{t}*x + tp*x);
        P= P+ tp;
        
             
      d = diag(P);
d1 = diag(P,1);
d2 = diag(P,n1);
d = diag(d);
d1 = diag(d1,1);
d2 = diag(d2,n1);
U2 = d +d1 +d2;
d1t = diag(P,-1);
d2t = diag(P,-n1);
ratios = 1./(diag(d));
d1t = ratios(1:N-1).*d1t;
d2t = ratios(1:N-n1).*d2t;
L2 = In+diag(d1t,-1)+diag(d2t,-n1);
% 
    invp = @(x) FFP(L2,U2,x);
        
        
        
%         invp = @(x) fastinvP(P,x);
        invP =  A_operator(@(x) invp(x), @(x) invp(x));
        
        [y(:,t),k1,res(:,t),err(:,t)] =  PCG_operator(S,Atb(:,t),invP,input.cgiter,y(:,t),1e-15,1,y(:,t));
    end
    y = project(y);

    output.err(:,itr) = sum(err,2);
    output.res(:,itr) = sum(res,2);
    output.snr(iterno)=snr(y, input.f);
    output.xtime(iterno)=cputime-t00;
%     output.funv(iterno)=get_f(A, y, alpha, b, f,n1,n2);
    
    
    if norm(yp-y)/norm(yp) < input.tol
        break
    end
end

output.res = mean(output.res,2);
output.err = mean(output.err,2);
output.y=y; %output.x=z;


function [f, f_prev]=get_f(A, y, alpha, b, f,n1,n2)
v1 = alpha*TV(y,n1,n2);
%         v2 = beta*norm(Phi*y,1);
r3 = A*y - b; v3 = 0.5*(norm(r3).^2);
f_prev = f;
f = v1 + v3;
return

function TVx=TV(x,n1,n2)
% Total Variation norm of x, x is a n by n matrix
TV_eps = 0;
grad_x = [Grad1(x,n1,n2) Grad2(x,n1,n2)];
if 1
    pt_sqsum = sum(grad_x.*grad_x,2);
    if TV_eps == 0; TVx = sum(sqrt(pt_sqsum)); else TVx = sum(sqrt(pt_sqsum+TV_eps)); end
else
    TVx = norm(grad_x(:),1);
end
return
function p=Grad1(u,n1,n2)
% backward finite difference along dim 1
u = reshape(u,n1,n2);
p = [u(1,:);diff(u,1,1);];
p = p(:);
return

function q=Grad2(u,n1,n2)
% backward finite difference along dim 2
u = reshape(u,n1,n2);
q = [u(:,1) diff(u,1,2)];
q = q(:);
return
