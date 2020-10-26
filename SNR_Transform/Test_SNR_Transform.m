clc;clear;
n=128*128; T=4;
ref=randn(n,T); 
est=ref+0.01*randn(n,T);
est1=ref+0.01*randn(n,T);

for t=1:T,
    S(1,t)=snr(est(:,t), ref(:,t));
    S(2,t)=snr(est1(:,t), ref(:,t));
end
S0(1,1)=snr(est(:), ref(:));
S0(2,1)=snr(est1(:), ref(:));

snr0=MT_SNR(S, ref);

a=find(abs(snr0-S0)>1e-6)
