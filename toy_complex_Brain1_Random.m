%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% This script is used to obtain the results in the paper
%%%% "Calibrationless Parallel MRI with Joint Total Variation
%%%% Regularization"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;
rand('state',0); rand('state',0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load rawdata_brain;
im = ifft2c(raw_data);

% Parameter Settings

sigma = 0.01;     % data noise
alpha = 4e-2; beta = 4e-2;  % regularization parameters
d = 4;   %parameter to control sampling ratio  d=3.45 for 33% d=4 for 25% d=4.4 for 20%
maxIter=50; %iteration number for CSSENSE, SPGL1 and JTVMRI

m=256; n=m;
coils = 8; T=coils;


%% Generate Mask
samplingratio = 0;
OMEGA = RandMask_rect(double(m/d),double(n/d),m,n);
mask = RandMask_InverseTransfer(OMEGA,m,n);
%mask(m/2-14:m/2+15,n/2-14:n/2+15) = ones(30,30);  % Why? The central white box? 
OMEGA = find(fftshift(mask)==1);
k = 2*length(OMEGA)+1;
samplingratio = samplingratio + length(OMEGA);
omask = zeros(m, n);
omask(OMEGA) = 1;

figure;
imshow(omask, []);
figure;
imshow(mask, []);