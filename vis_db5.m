clear all; close all;
db_pars = struct();
db_pars.resize = [256, 256];
db_pars.normalize = 1;
[data, fdata] = load_db(5, db_pars);
[m, n, T] = size(fdata);

im = sos(abs(data)); 
figure; imshow(im, []); colormap(gray(256));
im1 = abs(data(:,:,1));
im2 = abs(data(:,:,4));
im3 = abs(data(:,:,8));
fim1 = abs(fdata(:,:,1));
fim2 = abs(fdata(:,:,4));

% figure; imshow(im1, []); colormap(gray(256)); 
% figure; imshow(im2, []);
% figure; colormap((gray(256))); imshow(fim1, []);
% figure; colormap((gray(256))); imshow(fim2, []);
% 
im_tv1(:, :, 1) = conv2(im1, [-1 1], 'same');
im_tv1(:, :, 2) = conv2(im1, [-1 1]', 'same');
figure; imshow(sos(im_tv1), []);
% 
im_tv2(:, :, 1) = conv2(im2, [-1 1], 'same');
im_tv2(:, :, 2) = conv2(im2, [-1 1]', 'same');
figure; imshow(sos(im_tv2), []);
im_tv3(:, :, 1) = conv2(im3, [-1 1], 'same');
im_tv3(:, :, 2) = conv2(im3, [-1 1]', 'same');
figure; imshow(sos(im_tv3), []);
% 
% mask_pars = struct();
% mask_pars.image_size = [m, n];
% mask_pars.central_window = [30, 30];
% mask_pars.d = 4;
% mask = load_mask('random', mask_pars);
% figure;imshow(mask,[]);