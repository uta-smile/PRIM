function output = FIRLS_TV(input,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% minimize alpha*||x||_TV + 0.5*||Ax-b||_2^2

% Related papers:
%
% Chen Chen, Yeqing Li, Leon Axel and Junzhou Huang, "Real Time Dynamic MRI with Dynamic Total Variation", the Annual International Conference on Medical Image Computing and Computer Assisted Intervention (MICCAI), 2014.
%
% Chen Chen, Yeqing Li, Leon Axel and Junzhou Huang, "Real Time Dynamic MRI by Exploiting Spatial and Temporal Sparsity with Dynamic Total Variation", submitted to IEEE Transactions on Image Processing.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ratio =  input.ratio;
n1=input.n1; n2=input.n2; N=n1*n2;
A=input.A; b=input.b;
At = input.At;

Atb = At*b;

y= Atb;
if (isfield(input,'x0'))
    y = input.x0;
end

alpha = input.alpha;

if isnumeric(A)
    AtA=At*A;
else
    AtA =  A_operator(@(x) At*reshape(A*x,[n1,n2]), @(x) At*(A*x));
end
Q1 = input.Q1;
Q2 = input.Q2;
In = speye(N);
f0 = input.f;
for itr = 1:input.no    % total iter counter       
    % update weight
    D1 = zeros(size(y(:,1)));
    D2 = D1;
    
    D1=D1 + (abs(Q1*y)+eps).^(2);
    D2=D2 + (abs(Q2*y)+eps).^(2);
    
    D = (D1+D2).^(-0.5); %weight matrix
    
    P=ratio*In;
    tp=[(1:N)' (1:N)' D];
    dd=spconvert(tp);
    tp =  alpha*((Q1'*dd*Q1)+(Q2'*dd*Q2));
    P= P+ tp; %preconditioner
    S = A_operator(@(x) AtA*x + tp*x, @(x) AtA*x + tp*x); %system matrix
    
    
    %construct LU decompostion
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
    invP =  A_operator(@(x) invp(x), @(x) invp(x));
    
    % PCG
    [y,k1,res,err] =  PCG_operator(S,Atb,invP,input.cgiter,y,1e-15,1,f0(:)); 
end

output.y=y; 


