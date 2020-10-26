% This script performs a thresholding experiment as described in the paper
% Sparse MRI: The application of compressed sensing for rapid MR Imaging.
%
% A transform of an image is calculated, the coefficients below a threshold are
% zeroed, and an image is reconstructed from resulting subset of the largest transform coefficients


% load brain image

if exist('img')==0
	disp('loading brain image')
	load brain
	img = img./max(abs(img(:)));
end

[nx,ny] = size(img);
IMSIZE = nx*ny;

% compute DCT, wavelet and finite-difference transform
IMG_dct = FDCT(img,8);
IMG_wav = FWT2_PO(img,3,MakeONFilter('Symmlet',4));
IMG_d = D(img);

% sort coefficient from large to small
idx_dct = sort(abs(IMG_dct(:)),1,'descend');
idx_wav = sort(abs(IMG_wav(:)),1,'descend');
idx_d = sort(abs(IMG_d(:)),1,'descend');


c = 0;

% iterate on several threshold parameters and show results.
for pctg = floor([0.01,0.05,0.1,0.2,0.3,0.5]*IMSIZE);
c = c+1;

% threshold DCT coefficients
thresh = idx_dct(pctg);
tmp = IMG_dct(:);
tmp(find(abs(tmp)<thresh))=0;
rec_dct(:,:,c) = IDCT(tmp,8,nx,ny);

% threshold wavelet coefficients
thresh = idx_wav(pctg);
tmp = IMG_wav;
tmp(find(abs(tmp)<thresh))=0;
rec_wav(:,:,c) = IWT2_PO(tmp,3,MakeONFilter('Symmlet',4));

% threshold finite differences
thresh = idx_d(pctg);
tmp = IMG_d;
tmp(find(abs(tmp)<thresh))=0;
rec_d(:,:,c) = invD(tmp,[nx,ny]);


end

tmp = [];
for n=1:5
tmp = cat(2,tmp,abs(cat(1, rec_dct(:,:,n), rec_wav(:,:,n), rec_d(:,:,n))));
end
figure, imshow(tmp,[],'InitialMagnification',100), title('1% , 5% , 10% , 20% , 30% , 50% '), 
ylabel('        Finite Diff.                  Wavelet                     DCT '), drawnow,

disp('Done')



