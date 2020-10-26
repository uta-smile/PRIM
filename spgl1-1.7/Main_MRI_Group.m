clc;clear;
% Initialize random number generator
randn('state',0); 

load MT_Phantom; X=MT_Phantom(:,:,1:2);
m=128; n=m;
T=size(X, 3);
for t=1:T,
    temp=imresize(X(:,:,t), [m n]);
    temp=im256(temp);
    F(:,:,t)=temp;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Preparing the data and operators
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[m n, T] = size(F); N = m*n; f=F(:); pn=m*n;
d = 4;%dd(iter);   %%%% choose 4, 4.4, 5, 6.3

wav = daubcqf(2);W = @(x) midwt(x,wav);WT = @(x) mdwt(x,wav);
Phi= A_operator(@(x) WT(x), @(x) W(x));    % notice Phi=WT

for t=1:T,
     [OMEGA{t}] = RandMask_rect(double(m/d),double(n/d),m,n);
     k(t) = 2*length(OMEGA{t})+1;
end
for t=1:T,
     start=1+sum(k(1:t-1)); stop=sum(k(1:t));
     sz(t,:)=[start, stop, pn];
end

for t=1:T,
     f2=F(:,:, t); x0((t-1)*pn+1:t*pn,1)=WT(f2(:)); 
     f0((t-1)*pn+1:t*pn,1)=f2(:);
     R{t} = @(x) A_fhp_rect(x, OMEGA{t}, m, n);
     RT{t} = @(x) At_fhp_rect(x, OMEGA{t}, m, n);
     AO{t} = @(x) R{t}(W(x)); AOT{t} = @(x) WT(RT{t}(x)); 
     A{t} = A_operator(@(x) AO{t}(x), @(x) AOT{t}(x));
     sigma = 0.01; noise = sigma*randn(k(t),1); 
     start=sz(t,1); stop=sz(t,2);
     b(start:stop,1) = A{t}* x0((t-1)*pn+1:t*pn,1)+ noise;
end

nGroups = pn; 
temp=[1:pn]';temp=repmat(temp, [1, T]);
groups=temp(:);
    
% Solve unweighted version
opts = spgSetParms('verbosity',0);
opts.Phi=Phi;opts.f0=f0; opts.sz=sz; opts.iterations=50;
[x,r,g,info] = spg_mt(A,b,groups,0,opts);
mse=norm(x-x0,2)/norm(x0,2)
figure; 
plot(info.Trace_Time, info.Trace_SNR, 'b-', 'linewidth', 2)