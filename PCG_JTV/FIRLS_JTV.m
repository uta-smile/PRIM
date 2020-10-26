%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FIRLS for total variation (TV) minimization
%%% minimize_y 0.5*||Ay-b|| + lambda*||D y||_{2,1}

% %%% Input:
% input.A: the sensing matrix
% input.b: the measurement
% input.D1, input.D2: discrete finite difference matrices for TV
% input.lambda: regularization parameter
% input.cgiter: number of inner CG iterations
% input.f: the groudtruth. Used for calculating RMSE etc.
% input.ratio: the mean value of AtA
% input.tol: stopping tolerance
% input.l: the lower bound of x, may be 0 or -inf
% input.u; the upper bound of x, may be 1, 255 or inf
% input.n1, input.n2: the size of the image

% %%% Output
% output.y: the reconstructed image
% output.rel: the relative error
% output.snr: the SNR
% output.xtime: CPU time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Contact
%%%% Chen Chen (chenchen.cn87@gmail.com)
%%%% Junzhou Huang (jzhuang@uta.edu) University of Texas at Arlington

%%%% Related Papers
%%%% Chen Chen, Junzhou Huang, Lei He, and Hongsheng Li. "Preconditioning for accelerated iteratively reweighted least squares in structured sparsity reconstruction."
%%%% In IEEE Conference on Computer Vision and Pattern Recognition (CVPR) , pp. 2713-2720. IEEE, 2014.

%%%% Chen, Chen, Junzhou Huang, Lei He, and Hongsheng Li. "Fast Iteratively Reweighted Least Squares Algorithms for Analysis-Based Sparsity Reconstruction." arXiv preprint arXiv:1411.5057 (2014).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function output = FIRLS_JTV(input,varargin)

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

m=input.n1; n=input.n2; N=m*n;
A=input.A; b=input.b;
T = size(b,2);

for t=1:T
    Atb(:,t)=A{t}'*b{t};
end

y= Atb;
% y = zeros(size(Atb));
lambda = input.lambda;

if isnumeric(A{t})
    for t=1:T
        AtA{t}=A{t}'*A{t};
    end
else
    for t=1:T
        AtA{t} =  A_operator(@(x) A{t}'*(A{t}*x), @(x) A{t}'*(A{t}*x));
    end
end

D1 = input.D1;
D2 = input.D2;
In = speye(N);

output.funval(1) = CalculateJTVFunctionVal(A, reshape(y, [m, n, T]), b, lambda);
output.rmse(1)=RMSE(y, input.f);
output.snr(1)=snr(y(:), input.f(:));
t0 = tic;
eps = 1e-10;

for itr = 1:input.no    % total iter counter
    
    yp=y;
    Q1 = zeros(size(y(:,1)));
    Q2 = Q1;
    for t =1:T
        Q1=Q1 + (abs(D1*y(:,t))+eps).^(2);
        Q2=Q2 + (abs(D2*y(:,t))+eps).^(2);
    end
    Q = (Q1+Q2).^(-0.5);
    
    tp=[(1:N)' (1:N)' Q];
    W=spconvert(tp); %weight
    tp =  lambda*((D1'*W*D1)+(D2'*W*D2));
    
    for t=1:T
        P=ratio(t)*In;
        P= P+ tp;
        S = A_operator(@(x) AtA{t}*x + tp*x, @(x) AtA{t}*x + tp*x);
        
        [L,U] = inexactLUdecomposition(P, m, N);
        invp = @(x) FFP(L,U,x);
        invP =  A_operator(@(x) invp(x), @(x) invp(x));
        
        y(:,t) =  PCG_operator(S,Atb(:,t),invP,input.cgiter,y(:,t),1e-15,1,y(:,t));
    end
    y = project(y);
    
    output.xtime(itr)=toc(t0);
    
    output.rel(itr)=norm(y-yp, 'fro')/norm(yp, 'fro');
    output.snr(itr+1)=snr(y(:), input.f(:));
    output.rmse(itr+1)=RMSE(y, input.f);
    
    output.funval(itr+1) = CalculateJTVFunctionVal(A, reshape(y, [m, n, T]), b, lambda);
    
    
    if(output.rel(end)<input.tol)
        display(['Done at iteration ' num2str(itr)]);
        break;
    end
end

output.y=reshape(y, [m, n, T]);

function [L,U] = inexactLUdecomposition(P, n1, N)

d = diag(P);
d1 = diag(P,1);
d2 = diag(P,n1);
d = diag(d);
d1 = diag(d1,1);
d2 = diag(d2,n1);
U = d +d1 +d2;
d1t = diag(P,-1);
d2t = diag(P,-n1);
ratios = 1./(diag(d));
d1t = ratios(1:N-1).*d1t;
d2t = ratios(1:N-n1).*d2t;
In = speye(N);
L = In+diag(d1t,-1)+diag(d2t,-n1);

function val = CalculateJTVFunctionVal(A, X, b, lambda) 
    r = 0;
    [m,n, T]=size(X);
    PS1 = zeros(size(X));
    PS2 = zeros(size(X));
    for t = 1:T
        r = r + norm(abs(A{t}*reshape(X(:, :, t), [m*n, 1]) - b{t}), 'fro'); 
    end
    PS1(1:m-1, :, :) = X(1:m-1,:,:)-X(2:m,:,:);
    PS2(:, 1:n-1, :) = X(:,1:n-1,:)-X(:,2:n,:);
    tv = sqrt(abs(sum(PS1, 3)).^2+abs(sum(PS2, 3)).^2);
    val = r + sum(tv(:));

