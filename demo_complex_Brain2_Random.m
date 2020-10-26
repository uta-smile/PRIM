%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is an example for parrallel MRI 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Related papers:
% 
% Chen Chen, Yeqing Li, and Junzhou Huang, "Calibrationless Parallel MRI with Joint Total Variation Regularization", the Annual International Conference on Medical Image Computing and Computer Assisted Intervention (MICCAI), 2013.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;

load brain_8ch;
im_Full = ifft2c(DATA);

alpha = 1e-2; % regularization parameters
maxIter=400;

[m n, T] = size(im_Full); N = m*n;

% prepare sampling mask
numlines = 70; % change for different sampling ratios
OMEGA=fftshift(MRImask(m,numlines));  
[mask] = RandMask_InverseTransfer(OMEGA,m,n);
mask(m/2-14:m/2+15,n/2-14:n/2+15) = ones(30,30);
OMEGA = find(fftshift(mask)==1);
figure;imshow(mask,[]);

for t=1:T
    [masks{t}] = mask;
end

%generate Fourier sampling operator
for t=1:T
A{t} = p2DFT(masks{t}, size(im_Full(:,:,1)), 1, 2);
b{t} = A{t}*im_Full(:,:,t);
im_dc(:,:,t) = A{t}'*b{t};
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% JTV
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input.f=im_Full; input.n1=m;input.n2=n;
input.alpha=alpha;
input.L=1;
input.no=maxIter; %%%%  Change no for different iteration numbers

fprintf('calling the function JTV.....\n');
out = FISTA_JTV(b, A, input);
im_rec = out.y;

im_JTVMRI = sos(im_rec);
im_FullSoS = sos(im_Full);
im_dcSoS = sos(im_dc);