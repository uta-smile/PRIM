% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L+S reconstruction of dynamic contrast-enhanced abdominal MRI acquired
% with golden-angle radial sampling
% 
% Temporal resolution is flexible and determined by the user in the 
% variable nspokes 
%
% Ricardo Otazo (2013)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;close all;
addpath('nufft_toolbox/');
% number of spokes to be used per frame (Fibonacci number)
nspokes=21;
% load radial data
load abdomen_dce_ga.mat
[nx,ny,nc]=size(b1);
[nr,ntviews,nc]=size(kdata);
% number of frames
nt=floor(ntviews/nspokes);
% crop the data according to the number of spokes per frame
kdata=kdata(:,1:nt*nspokes,:);
k=k(:,1:nt*nspokes);
w=w(:,1:nt*nspokes);
% sort the data into a time-series of undersampled images
for ii=1:nt
    kdatau(:,:,:,ii)=kdata(:,(ii-1)*nspokes+1:ii*nspokes,:);
    ku(:,:,ii)=k(:,(ii-1)*nspokes+1:ii*nspokes);
    wu(:,:,ii)=w(:,(ii-1)*nspokes+1:ii*nspokes);
end
% multicoil NUFFT operator
param.E=MCNUFFT(ku,wu,b1);
param.d=kdatau;
recon_nufft=param.E'*param.d;
clear kdata k ku wu w
% L+S reconstruction ******************************************************
param.lambda_L=0.025;
param.lambda_S=0.5*max(abs(recon_nufft(:)));
param.nite=20;
param.tol=0.0025;
[L,S] = lps_tv(param);
L=flipdim(L,1);
S=flipdim(S,1);
LplusS=L+S;

% display 4 frames
LplusSd=LplusS(65:336,:,1);LplusSd=cat(2,LplusSd,LplusS(65:336,:,9));LplusSd=cat(2,LplusSd,LplusS(65:336,:,16));LplusSd=cat(2,LplusSd,LplusS(65:336,:,25));
Ld=L(65:336,:,1);Ld=cat(2,Ld,L(65:336,:,9));Ld=cat(2,Ld,L(65:336,:,16));Ld=cat(2,Ld,L(65:336,:,25));
Sd=S(65:336,:,1);Sd=cat(2,Sd,S(65:336,:,9));Sd=cat(2,Sd,S(65:336,:,16));Sd=cat(2,Sd,S(65:336,:,25));

figure;
subplot(3,1,1),imshow(abs(LplusSd),[0,5e-4]);ylabel('L+S')
subplot(3,1,2),imshow(abs(Ld),[0,5e-4]);ylabel('L')
subplot(3,1,3),imshow(abs(Sd),[0,5e-4]);ylabel('S')








