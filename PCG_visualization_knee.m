close all; clear;

load PCG_JTV_db5_fused.mat

% cx = 60; cy = 120;
% cx = 90; cy = 170;
% cx = 130; cy = 100;
cx = 135; cy = 130;

ims(:,:,1) = abs(sos(data));
ims(:,:,2) = abs(sos(im_grappa));
ims(:,:,3) = abs(sos(im_cgSPIRiT));
ims(:,:,4) = abs(sos(im_cssense));
ims(:,:,5) = abs(sos(im_sake));
ims(:,:,6) = abs(sos(im_fista_jtv));
ims(:,:,7) = abs(sos(im_pcg_jtv));

for i = 1:7
    ims(:, :, i) = im256(ims(:, :, i));
end

zoomedim = zoominim(ims,1,cx,cy);

subplot(2,4,1); imshow(zoomedim(:,:,1),[]); title('Original');
subplot(2,4,2); imshow(zoomedim(:,:,2),[]); title('GRAPPA');
subplot(2,4,3); imshow(zoomedim(:,:,3),[]); title('CGSPIRiT');
subplot(2,4,4); imshow(zoomedim(:,:,4),[]); title('CSSENSE');
subplot(2,4,5); imshow(zoomedim(:,:,5),[]); title('SAKE');
subplot(2,4,6); imshow(zoomedim(:,:,6),[]); title('FISTA\_JTV');
subplot(2,4,7); imshow(zoomedim(:,:,6),[]); title('Proposed');

figure; imshow(zoomedim(:, :, 1), []); 
figure; imshow(zoomedim(:,:,2),[]); 
figure; imshow(zoomedim(:,:,3),[]); 
figure; imshow(zoomedim(:,:,4),[]); 
figure; imshow(zoomedim(:,:,5),[]); 
figure; imshow(zoomedim(:,:,6),[]); 
figure; imshow(zoomedim(:,:,6),[]); 