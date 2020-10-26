function [X_den,P1,P2,iter]=denoise_TV_MT(Xobs,lambda,l,u,P1_init, P2_init,pars)
%This function implements the FISTA method for JTV denoising problems. 
%
% INPUT
% Xobs ..............................observed noisy images. [H, W, T]
% lambda ........................ parameter
% pars.................................parameters structure
% pars.MAXITER ..................... maximum number of iterations
%                                                      (Default=100)
% pars.epsilon ..................... tolerance for relative error used in
%                                                       the stopping criteria (Default=1e-4)
% pars.print ..........................  1 if a report on the iterations is
%                                                       given, 0 if the  report is silenced
% pars.tv .................................. type of total variation
%                                                      penatly.  'iso' for isotropic (default)
%                                                      and 'l1' for nonisotropic
%  
% OUTPUT
% X_den ........................... The solution of the problem 
%                                            min{||X-Xobs||^2+2*lambda*TV(X)}


%Define the Projection onto the box
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

% Assigning parameres according to pars and/or default values
flag=exist('pars', 'var');
if (flag&&isfield(pars,'MAXITER'))
    MAXITER=pars.MAXITER;
else
    MAXITER=100;
end
% if (flag&&isfield(pars,'epsilon'))
%     epsilon=pars.epsilon;
% else
%     epsilon=1e-4;
% end
% if(flag&&isfield(pars,'print'))
%     prnt=pars.print;
% else
%     prnt=1;
% end
if(flag&&isfield(pars,'tv'))
    tv=pars.tv;
else
    tv='iso';
end

[m,n,T]=size(Xobs);
% clear P; clear R;
if(isempty(P1_init))
    for t=1:T
        P1{t}=zeros(m-1,n);    P2{t}=zeros(m,n-1);
        R1{t}=zeros(m-1,n);    R2{t}=zeros(m,n-1);
    end
else
    for t=1:T
        P1{t}=P1_init{t};    P2{t}=P2_init{t};
        R1{t}=P1_init{t};    R2{t}=P2_init{t};
    end
end
tk=1;tkp1=1;count=0;i=0;

D=gpuArray.zeros(m,n,T);%fval=inf;fun_all=[];
while((i<MAXITER)&&(count<5))
%    fold=fval;  
    i=i+1;    
%     Dold=D;    
    Pold1=P1;Pold2=P2;    
    tk=tkp1;
    
    for t=1:T
        temp=Xobs(:,:,t)-lambda*Lforward(R1{t}, R2{t}, m, n);
        D(:,:,t)=project(temp);        
    end
    
    %%%%%%%%%%
    % Taking a step towards minus of the gradient
    step0=1/(8*lambda);
    for t=1:T
        [tQ1, tQ2]=Ltrans(D(:,:,t), m, n);
        P1{t}=R1{t}+step0*tQ1{1};
        P2{t}=R2{t}+step0*tQ2{1};
    end
    
    %%%%%%%%%%
    % Peforming the projection step
    switch tv
        case 'iso'
            A=[abs(P1{1}).^2;gpuArray.zeros(1,n)]+[abs(P2{1}).^2,gpuArray.zeros(m,1)];
            for t=2:T,
                A=A+[abs(P1{t}).^2;gpuArray.zeros(1,n)]+[abs(P2{t}).^2,gpuArray.zeros(m,1)];
            end
            A=sqrt(max(A,1));
            for t=1:T,
                P1{t}=P1{t}./A(1:m-1,:); 
                P2{t}=P2{t}./A(:,1:n-1);
            end
        otherwise
            error('unknown type of total variation. should be iso');
    end

    %%%%%%%%%%
    %Updating R and t
    tkp1=(1+sqrt(1+4*tk^2))/2;
    
    step=(tk-1)/tkp1;
    for t=1:T,
        R1{t}=P1{t}+step*(P1{t}-Pold1{t});
        R2{t}=P2{t}+step*(P2{t}-Pold2{t});
    end
    
%     re=norm(D-Dold,'fro')/norm(D,'fro');
%     if (re<epsilon)
%         count=count+1;
%     else
%         count=0;
%     end
%     C=Xobs-lambda*Lforward(P, m, n, T);
%     PC=project(C);
%     fval=-norm(C-PC,'fro')^2+norm(C,'fro')^2;
%     fun_all=[fun_all;fval];
end
X_den=D;iter=i;


function X=Lforward(P1, P2, m, n)

      X=gpuArray.zeros(m,n);
      X(1:m-1,:)=P1;
      X(:,1:n-1)=X(:,1:n-1)+P2;
      X(2:m,:)=X(2:m,:)-P1;
      X(:,2:n)=X(:,2:n)-P2;
   end
 
   function [P1, P2]=Ltrans(X, m, n)
%       [m,n]=size(X);
      
          P1{1}=X(1:m-1,:,:)-X(2:m,:,:);
          P2{1}=X(:,1:n-1,:)-X(:,2:n,:);
      
   end
end