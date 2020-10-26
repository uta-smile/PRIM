function [ res, fres ] = runPocsSPIRiT( fdata, mask, pars )
%RUNCGSPIRIT Summary of this function goes here
%   Detailed explanation goes here


kSize = pars.kSize;  % SPIRiT kernel size
nIter = pars.nIter; % number of iteration; phantom requires twice as much as the brain.
CalibTyk = pars.CalibTyk;  % Tykhonov regularization in the calibration
wavWeight = pars.wavWeight;  % Wavelet soft-thresholding regularization in the reconstruction (SPIRiT only)

[CalibSize, dcomp] = getCalibSize(mask);  % get size of calibration area from mask
[m, n, T] = size(fdata);
fdata = fdata.*repmat(mask, [1,1,T]); % multiply with sampling matrix

% DATAcomp = fdata.*repmat(dcomp,[1,1,T]);
% scale_fctr = norm(DATAcomp(:))/sqrt(T)/20;
% fdata = fdata/scale_fctr;
% DATAcomp = DATAcomp/scale_fctr;

% im_dc = ifft2c(DATAcomp);
% im = ifft2c(fdata); 
% im = im/scale_fctr;

disp('performing calibration for SPIRiT')
kCalib = crop(fdata,[CalibSize,T]);
kernel = zeros([kSize,T,T]);

[AtA] = corrMatrix(kCalib,kSize);
for t=1:T
	kernel(:,:,:,t) = calibrate(AtA,kSize,T,t,CalibTyk);
end
GOP = SPIRiT(kernel, 'fft',[m ,n]);

disp('performing pocs reconstruction')
tic;
% [fres] = pocsSPIRiT(fdata,GOP,nIter, fdata,wavWeight,0);
[fres] = pocsSPIRiT(fdata,GOP,nIter, zeros(size(fdata)),wavWeight,0);
toc

res = ifft2c(fres);


end

