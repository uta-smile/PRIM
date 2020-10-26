function snr0=MT_SNR(SNR, Ref);
%%%% Obtain the global SNR from the SNRs of multiple components
%%%% Ref p x T where T is the task number and p is pixel number
%%%% SNR k x T, where k is iteration number and T is task number
[K, T]=size(SNR);
if size(Ref, 2)~=T
    return
end
for t=1:T
    v2(1, t)=var(Ref(:,t), 1);
end
v0=var(Ref(:),1);

for k=1:K,
    temp=v2./(10.^(SNR(k,:)/10));
   temp1=mean(temp);
   snr0(k,:)=10*log10(v0/temp1);
end
