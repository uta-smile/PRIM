% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L+S reconstruction of undersampled multicoil cardiac perfusion MRI
%
% Ricardo Otazo (2013)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;close all;
% load undersampled data 
load cardiac_perf_R8.mat;
[nx,ny,nt,nc]=size(kdata);
% L+S reconstruction ******************************************************
param.E=Emat_xyt(kdata(:,:,:,1)~=0,b1);
param.d=kdata;
param.T=TempFFT(3);
param.lambda_L=0.01;
param.lambda_S=0.01;
param.nite=50;
param.tol=0.0025;
[L,S] = lps_ist(param);
LplusS=L+S;

% display 4 frames
LplusSd=LplusS(33:96,33:96,2);LplusSd=cat(2,LplusSd,LplusS(33:96,33:96,8));LplusSd=cat(2,LplusSd,LplusS(33:96,33:96,14));LplusSd=cat(2,LplusSd,LplusS(33:96,33:96,24));
Ld=L(33:96,33:96,2);Ld=cat(2,Ld,L(33:96,33:96,8));Ld=cat(2,Ld,L(33:96,33:96,14));Ld=cat(2,Ld,L(33:96,33:96,24));
Sd=S(33:96,33:96,2);Sd=cat(2,Sd,S(33:96,33:96,8));Sd=cat(2,Sd,S(33:96,33:96,14));Sd=cat(2,Sd,S(33:96,33:96,24));

figure;
subplot(3,1,1),imshow(abs(LplusSd),[0,1]);ylabel('L+S')
subplot(3,1,2),imshow(abs(Ld),[0,1]);ylabel('L')
subplot(3,1,3),imshow(abs(Sd),[0,1]);ylabel('S')