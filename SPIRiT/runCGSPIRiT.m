function [ res, fres ] = runCGSPIRiT( fdata, mask, pars )
%RUNCGSPIRIT Summary of this function goes here
%   Detailed explanation goes here


kSize = pars.kSize;  % SPIRiT kernel size
nIterCG = pars.nIterCG; % number of iteration; phantom requires twice as much as the brain.
% mask_type = 'unif'; % options are: 'unif','random4','random3'
CalibTyk = pars.CalibTyk;  % Tykhonov regularization in the calibration
ReconTyk = pars.ReconTyk;  % Tykhovon regularization in the reconstruction (SPIRiT only)

[m, n, T] = size(fdata);
% im = ifft2c(fdata);

[CalibSize, dcomp] = getCalibSize(mask);
DATA = fdata.*repmat(mask,[1, 1, T]); % multiply with sampling matrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% scale the data such that the zero-filled density compensated      %%%%%%%%%
% k-space norm is 1. This is useful in order to use similar         %%%%%%%%%
% regularization penalty values for different problems.             %%%%%%%%%

% DATAcomp = DATA.*repmat(dcomp,[1, 1, T]);
% scale_fctr = norm(DATAcomp(:))/sqrt(coils)/20;
% DATA = DATA/scale_fctr;
% DATAcomp = DATAcomp/scale_fctr;

% im_dc = ifft2c(DATAcomp);
% im = im/scale_fctr;


disp('performing calibration for SPIRiT')
kCalib = crop(DATA ,[CalibSize,T]);

kernel = calibSPIRiT(kCalib, kSize, T, CalibTyk);
GOP = SPIRiT(kernel, 'fft',[m, n]);

disp('performing CG reconstruction')
tic;
[fres, ~] = cgSPIRiT(DATA, GOP, nIterCG, ReconTyk, DATA);
toc

res = ifft2c(fres);


end

