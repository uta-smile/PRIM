function output = FISTA_JTV(b, A, input)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% minimize alpha*JTV(X) + beta*(||Phi*X||_{2,1}+ ||Phi*X||_{Forest}) + 0.5*sum_{s}||R_{s}X(:,s)-b_{s}||_2^2
%%% Phi:DWT,  Phi': IDWT
%%% Jan. 16, 2013, Written by Chen Chen at University of Texas at
%%% Arlington 

%%% Chen Chen, Yeqing Li, and Junzhou Huang, "Calibrationless Parallel MRI with Joint Total Variation Regularization", the Annual International Conference on Medical Image Computing and Computer Assisted Intervention (MICCAI), 2013.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parsin.MAXITER=1; parsin.tv='iso'; 

if iscell(b)
    T=length(b); 
else
    T = 1;
    A = {A};
    b = {b};
end

Lx=input.L; 
n1=input.n1; n2=input.n2; N=n1*n2;
alpha = input.alpha; 


for t=1:T,
    Atb{t}=A{t}'*b{t};
end

y=zeros(input.n1,input.n2,T);
% set initial value
for t = 1:T
    y(:,:,t) = Atb{t};
end
yr=zeros(size(y));
tnew=1;

output.funval(1) = CalculateJTVFunctionVal(A, y, b, alpha);
output.rmse(1)=RMSE(y, input.f);
output.snr(1)=snr(y(:), input.f(:));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t0 = tic;
for itr = 1:input.no  % total iter counter        
    told=tnew;
    yp=y;
       
    %update yg    
    for t=1:T
        temp=A{t}'*(A{t}*yr(:,:,t));
        temp=temp-Atb{t};
        yg(:,:,t)=yr(:,:,t)-temp/Lx; 
    end
    
    if alpha 
        if (itr==1)
            [y, P1, P2]=denoise_TV_MT(yg, alpha/Lx,-inf,inf,[],[], parsin);
        else
            [y, P1, P2]=denoise_TV_MT(yg, alpha/Lx,-inf,inf,P1, P2,parsin);
        end
    else
        y = yg;
    end
    
    tnew=(1+sqrt(1+4*told^2))/2; 
    yr=y+((told-1)/tnew)*(y-yp);  
    
    output.xtime(itr)=toc(t0); 
    
    output.snr(itr+1)=snr(y(:), input.f(:));
    output.funval(itr+1) = CalculateJTVFunctionVal(A, y, b, alpha);
    
    output.rmse(itr+1)=RMSE(y, input.f);
    
    
end

output.y=y;

function val = CalculateJTVFunctionVal(A, X, b, lambda) 
    r = 0; 
    [m,n, T]=size(X);
    PS1 = zeros(size(X));
    PS2 = zeros(size(X));
    for t = 1:T
        r = r + norm(abs(A{t}*X(:, :, t) - b{t}), 'fro'); 
    end
    PS1(1:m-1, :, :) = X(1:m-1,:,:)-X(2:m,:,:);
    PS2(:, 1:n-1, :) = X(:,1:n-1,:)-X(:,2:n,:);
    tv = sqrt(abs(sum(PS1, 3)).^2+abs(sum(PS2, 3)).^2);
    val = r + sum(tv(:));
    