function output = FISTA_pJTV_GPU(b, A, pars)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% minimize alpha*JTV(X) + beta*(||Phi*X||_{2,1}+ ||Phi*X||_{Forest}) + 0.5*sum_{s}||R_{s}X(:,s)-b_{s}||_2^2
%%% Phi:DWT,  Phi': IDWT
%%% Jan. 16, 2013, Written by Chen Chen at University of Texas at
%%% Arlington 

%%% Chen Chen, Yeqing Li, and Junzhou Huang, "Calibrationless Parallel MRI with Joint Total Variation Regularization", the Annual International Conference on Medical Image Computing and Computer Assisted Intervention (MICCAI), 2013.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parsin.MAXITER=20; parsin.tv='iso'; 

if iscell(b)
    T=length(b); 
else
    T = 1;
    A = {A};
    b = {b};
end

Lx=pars.L; 
n1=pars.n1; n2=pars.n2; N=n1*n2;
lambda = pars.lambda; 


for t=1:T,
    Atb{t}=A{t}'*b{t};
end

y=zeros(pars.n1,pars.n2,T);
yr=zeros(size(y));
tnew=1;
t00 = cputime();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for itr = 1:pars.no  % total iter counter        
    told=tnew;
    yp=y;
       
    %update yg    
    for t=1:T
        temp=A{t}'*(A{t}*yr(:,:,t));
        temp=temp-Atb{t};
        yg(:,:,t)=yr(:,:,t)-temp/Lx; 
    end
    
%     if (itr==1)
%         [y, P1, P2]=denoise_TV_MT_GPU(yg, alpha/Lx,-inf,inf,[],[], parsin);
%     else
%         [y, P1, P2]=denoise_TV_MT_GPU(yg, alpha/Lx,-inf,inf,P1, P2,parsin);
%     end

    parsin.lambda = lambda;
    parsin.gamma  = pars.gamma;
    parsin.max_iter = pars.inner_max_iter;
    
    y = pJTV_GPU_complex(yg, parsin);
    
    tnew=(1+sqrt(1+4*told^2))/2; 
    yr=y+((told-1)/tnew)*(y-yp);    
    
    output.snr(itr)=snr(y(:), pars.f(:));

    output.xtime(itr)=cputime-t00; 

    output.rmse(itr)=RMSE(y, pars.f);
end

output.y=y;
end