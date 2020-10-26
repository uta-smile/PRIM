clear all; close all;
db_pars.resize = [128, 128];
% db_pars.normalize = 1;
[data, fdata] = load_db(2);
% regularization parameters
% alpha = 4e3; % for db 1
alpha = 1e-4; 

maxIter=100;
[m n, T] = size(fdata); N = m*n;

%% Random Mask Generation
mask_pars.image_size = [m, n];
mask_pars.central_window = [30, 30];
mask_pars.d = 4;
mask_pars.line_num = 70;
mask = load_mask('random', mask_pars);



%% JTV
A = {}; b = {};
for t=1:T
A{t} = p2DFT(mask, size(data(:,:,1)), 1, 2);
b{t} = A{t}*data(:,:,t);
im_dc(:,:,t) = A{t}'*b{t};
end
input.f=data; input.n1=m;input.n2=n;
input.alpha=alpha;
input.L=1;
input.no=maxIter; %%%%  Change no for different iteration numbers

fprintf('calling the function JTV.....\n');
t0 = tic;
out = FISTA_JTV(b, A, input); t1 = toc(t0);
im_jtv = out.y;


fprintf('RMSE: %.5f, time=%.2fs\n', RMSE(im_jtv, data), t1);

%% pJTV - OMP
A = {}; b = {};
for t=1:T
A{t} = p2DFT(mask, size(data(:,:,1)), 1, 2);
b{t} = A{t}*data(:,:,t);
im_dc(:,:,t) = A{t}'*b{t};
end
input.f=data; input.n1=m;input.n2=n;
input.lambda=alpha;
input.L=1;
input.no=maxIter; %%%%  Change no for different iteration numbers
input.inner_max_iter = 20;
input.gamma = 10;

fprintf('calling the function pJTV_OMP.....\n');
t0 = tic;
out = FISTA_pJTV_OMP(b, A, input); t2 = toc(t0);
im_pjtv_omp = out.y;


fprintf('RMSE: %.5f, time=%.2fs\n', RMSE(im_pjtv_omp, data), t2);

%% FISTA_JTV_GPU
% A = {}; b = {}; clear im_dc;
% gdata = gpuArray(data);
% gfdata = gpuArray(fdata);
% for t=1:T
% A{t} = p2DFT(mask, size(data(:,:,1)), 1, 2);
% b{t} = A{t}*gdata(:,:,t);
% im_dc(:,:,t) = A{t}'*b{t};
% end
% input.f=gdata; input.n1=m;input.n2=n;
% input.alpha=alpha;
% input.L=1;
% input.no=maxIter; %%%%  Change no for different iteration numbers
% 
% fprintf('calling the function JTV.....\n');
% t0 = tic;
% out = FISTA_JTV_GPU(b, A, input); t2 = toc(t0);
% im_jtv_gpu = out.y;
% fprintf('RMSE: %.5f, time=%.2fs\n', RMSE(im_jtv_gpu, gdata), t2);


%% FISTA_pJTV_GPU
% A = {}; b = {}; clear im_dc;
% gdata = gpuArray(data);
% gfdata = gpuArray(fdata);
% for t=1:T
% A{t} = p2DFT(mask, size(data(:,:,1)), 1, 2);
% b{t} = A{t}*gdata(:,:,t);
% im_dc(:,:,t) = A{t}'*b{t};
% end
% input.f=gdata; input.n1=m;input.n2=n;
% input.lambda=alpha;
% input.L=1;
% input.no=maxIter; %%%%  Change no for different iteration numbers
% input.inner_max_iter = 20;
% input.gamma = 10;
% 
% fprintf('calling the function JTV.....\n');
% t0 = tic;
% out = FISTA_pJTV_GPU(b, A, input); t3 = toc(t0);
% im_pjtv_gpu = out.y;
% fprintf('RMSE: %.5f, time=%.2fs\n', RMSE(im_pjtv_gpu, gdata), t3);


figure;imshow(abs(sos(im_jtv)),[]);
figure;imshow(abs(sos(im_pjtv_omp)),[]);
% 
% figure;imshow(mask,[]);
% figure; hold on;
% subplot(1,3,1); imshow(abs(sos(data)),[]); title('Ground Truth');
% subplot(1,3,2); imshow(abs(sos(im_jtv)),[]); title('JTVMRI');
% subplot(1,3,3  ); imshow(abs(sos(im_pjtv_omp)),[]); title('JTVMRI\_GPU');