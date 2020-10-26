function [res, kres] = runSAKE_ESPIRiT(b, pars)

ncalib = pars.ncalib;
ksize = pars.ksize; % ESPIRiT kernel-window-size

sakeIter = pars.sakeIter;
wnthresh = pars.wnthresh; % Window-normalized number of singular values to threshold
eigThresh_im = pars.eigThresh_im; % threshold of eigenvectors in image space
[sx,sy,Nc] = size(b);

calibc = crop(b,[ncalib,ncalib,8]);

%% Perform SAKE reconstruction to recover the calibration area

disp('Performing SAKE recovery of calibration');

tic;, calib = SAKE(calibc, [ksize], wnthresh,sakeIter, 0);toc

disp('Done')


%% Singular values of the calibration matrix and ESPIRiT Maps after SAKE
% Sake now shows a null space and improved Maps. 

[k,S] = dat2Kernel(calib,ksize);
[M,W] = kernelEig(k(:,:,:,1:floor(wnthresh*prod(ksize))),[sx,sy]);


%% Compute Soft-SENSE ESPIRiT Maps 
% crop sensitivity maps according to eigenvalues==1. Note that we have to
% use 2 sets of maps. Here we weight the 2 maps with the eigen-values

maps = M(:,:,:,end-1:end);

% Weight the eigenvectors with soft-senses eigen-values
weights = W(:,:,end-1:end) ;
weights = (weights - eigThresh_im)./(1-eigThresh_im).* (W(:,:,end-1:end) > eigThresh_im);
weights = -cos(pi*weights)/2 + 1/2;

% create and ESPIRiT operator
ESP = ESPIRiT(maps,weights);
nIterCG = 15; 


%% Reconsturctions
% ESPIRiT CG reconstruction with soft-sense and 2 sets of maps
disp('Performing ESPIRiT reconstruction from 2 maps')
tic; [kres, resESPIRiT] = cgESPIRiT(b, ESP, nIterCG, 0.01, b*0); toc
disp('Done');

res = ifft2c(kres);


end