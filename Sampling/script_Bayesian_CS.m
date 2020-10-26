%% Script to run Bayesian CS on complex-valued images

clear all;close all;
% addpath('utils')
% 
% load m2d_optsym_SL    % load undersampling masks, or generate new ones using the following cell
% 
% load SL_complex
F = double(imread('../pics/Heart.png'));F=imresize(F, [256 256]);
img=F;
im_size = size(img,1);
 L = size(img,3);
 
 
%  for t = 1:L
%     % undersampled k-space of imag(input)
%     y_imag(:,:,t) = -.5i*( fft2(img(:,:,t)).*m2d_sym(:,:,t) - circshift(flipud(fliplr(conj(fft2(img(:,:,t)).*m2d_sym(:,:,t)))),[1,1]) ) / im_size;
%     % undersampled k-space of real(input)
%     y_real(:,:,t) = .5*( fft2(img(:,:,t)).*m2d_sym(:,:,t) + circshift(flipud(fliplr(conj(fft2(img(:,:,t)).*m2d_sym(:,:,t)))),[1,1]) ) / im_size;
% end

%% generate symmetric sampling masks
% skip this cell if you want to use the previously generated mask stored in
% m2d_optsym_SL, otherwise run this cell to generate new sampling patterns

R = 5;   % acceleration factor
pdf = genPDF([1,256], 6, 1/R, 2, 0.1, 0);
for t = 1:L
    m2d_sym(:,:,t) = fftshift(repmat(genSampling(pdf,10,1),[im_size,1]));
end

figure;imshow(m2d_sym);
r = find(m2d_sym(:,:)==1); 
r =length(r);
r=r/256;
%% undersample in k-space

[k2,k1] = meshgrid(0:im_size-1,0:im_size-1);

fdx = ones(im_size,im_size) - exp(-2*pi*1i*k1/im_size);
fdy = ones(im_size,im_size) - exp(-2*pi*1i*k2/im_size);

for h = 1:L
    zf_real(:,:,h) = ifft2(y_real(:,:,h)) * im_size;
    zf_imag(:,:,h) = ifft2(y_imag(:,:,h)) * im_size;
        
    loc{2*h-1} = find(y_real(:,:,h));
    loc{2*h} = find(y_imag(:,:,h));
    
    y_redx = y_real(:,:,h).*fdx;
    y_redy = y_real(:,:,h).*fdy;
    y_imdx = y_imag(:,:,h).*fdx;
    y_imdy = y_imag(:,:,h).*fdy;
    
    tx{2*h-1} = cat(1, real(y_redx(loc{2*h-1})), imag(y_redx(loc{2*h-1})) );
    tx{2*h} = cat(1, real(y_imdx(loc{2*h})), imag(y_imdx(loc{2*h})) );
    
    ty{2*h-1} = cat(1, real(y_redy(loc{2*h-1})), imag(y_redy(loc{2*h-1})) );
    ty{2*h} = cat(1, real(y_imdy(loc{2*h})), imag(y_imdy(loc{2*h})) );

    img_x(:,:,2*h-1) = ifft2(fft2(real(img(:,:,h))).*fdx);
    img_x(:,:,2*h) = ifft2(fft2(imag(img(:,:,h))).*fdx);
    
    img_y(:,:,2*h-1) = ifft2(fft2(real(img(:,:,h))).*fdy);    
    img_y(:,:,2*h) = ifft2(fft2(imag(img(:,:,h))).*fdy);    
end

zf = zf_real+1i.*zf_imag;
tile( abs(zf),1,L,1), title(['Zf RMSE: ', num2str(100*norm(zf(:)-img(:))/norm(img(:))), ' percent'])


%% non-joint reconstruction, real and imag parts are reconstructed jointly,
%% but different constrasts are reconstructed separately

% this is for non-joint recon, the next cell does joint reconstruction

max_iter = 1e4;    % maximum number of iterations
eta = 1e-8;        % parameter for termination criterion
% It is best to run the iterations until the displayed error does not
% decrease any more. For this phantom, eta=1e-8 is too small, but it works
% well for in vivo data. max_iter is set to high value so that it
% terminates even if termination criterion is not met.

tic
for t = 1:L
     dx(:,2*t-1:2*t) = mt_CSfft2(im_size, {tx{2*t-1:2*t}}, {loc{2*t-1:2*t}}, 0, 0, eta, img_x(:,:,2*t-1:2*t), max_iter); 
     dy(:,2*t-1:2*t) = mt_CSfft2(im_size, {ty{2*t-1:2*t}}, {loc{2*t-1:2*t}}, 0, 0, eta, img_y(:,:,2*t-1:2*t), max_iter); 
end
toc

Dx = reshape( dx, [im_size,im_size,2*L] );  save Dx.mat Dx
Dy = reshape( dy, [im_size,im_size,2*L] );  save Dy.mat Dy


%% joint reconstruction, real and imag parts and different constrasts are
%% reconstructed jointly

max_iter = 1e4;
eta = 1e-8;

tic
     dx = mt_CSfft2(im_size, tx, loc, 0, 0, eta, img_x, max_iter); save dx_4task.mat dx;
     dy = mt_CSfft2(im_size, ty, loc, 0, 0, eta, img_y, max_iter); save dy_4task.mat dy;
toc

Dx = reshape( dx, [im_size,im_size,2*L] );
Dy = reshape( dy, [im_size,im_size,2*L] );
 

%% LS reconstruction from the gradient estimates
% after computing the spatial gradients with CS, find the image that
% matches these gradients and also the sampled k-space data


for h = 1:L
    Pls(:,:,h) = L2_image_from_edges(fft2(img(:,:,h)).*m2d_sym(:,:,h), Dx(:,:,2*h-1)+1i.*Dx(:,:,2*h), Dy(:,:,2*h-1)+1i.*Dy(:,:,2*h), 0);
end

tile([abs(Pls)],1,L,1), title(['Bayesian CS RMSE: ', num2str(100*norm(Pls(:)-img(:))/norm(img(:))), ' percent'])
tile([abs(img)],1,L,2)



