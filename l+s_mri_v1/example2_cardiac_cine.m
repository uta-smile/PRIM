% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L+S reconstruction of undersampled multicoil cardiac cine MRI
%
% Ricardo Otazo (2013)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;close all;
% load undersampled data 
load cardiac_cine_R6.mat;
[nx,ny,nt,nc]=size(kdata);
% L+S reconstruction ******************************************************
param.E=Emat_xyt(kdata(:,:,:,1)~=0,b1);
param.d=kdata;
param.T=TempFFT(3);
% correlated motion in the background
% param.lambda_L=0.0025;param.lambda_S=0.00125;
% stationary background
param.lambda_L=0.01;param.lambda_S=0.0025;
param.nite=50;
param.tol=0.0025;
[L,S] = lps_ist(param);
LplusS=L+S;

% display 4 frames
LplusSd=LplusS(65:192,65:192,2);LplusSd=cat(2,LplusSd,LplusS(65:192,65:192,8));LplusSd=cat(2,LplusSd,LplusS(65:192,65:192,14));LplusSd=cat(2,LplusSd,LplusS(65:192,65:192,20));
Ld=L(65:192,65:192,2);Ld=cat(2,Ld,L(65:192,65:192,8));Ld=cat(2,Ld,L(65:192,65:192,14));Ld=cat(2,Ld,L(65:192,65:192,20));
Sd=S(65:192,65:192,2);Sd=cat(2,Sd,S(65:192,65:192,8));Sd=cat(2,Sd,S(65:192,65:192,14));Sd=cat(2,Sd,S(65:192,65:192,20));

figure;
subplot(3,1,1),imshow(abs(LplusSd),[0,1]);ylabel('L+S')
subplot(3,1,2),imshow(abs(Ld),[0,1]);ylabel('L')
subplot(3,1,3),imshow(abs(Sd),[0,1]);ylabel('S')