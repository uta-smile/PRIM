%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is an example for parrallel MRI 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Related papers:
% 
% Chen Chen, Yeqing Li, and Junzhou Huang, "Calibrationless Parallel MRI with Joint Total Variation Regularization", the 
% Annual International Conference on Medical Image Computing and Computer Assisted Intervention (MICCAI), 2013.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;

db_pars.normalize = 0;
% db_pars.resize = [256, 256];
[data, fdata] = load_db(2, db_pars);

% regularization parameters
lambda_1 = 1; % JTV 
lambda_2 = 0; % Tensor-SVD 
max_iter = 100;

[m n, T] = size(fdata); N = m*n;

% Random Mask
mask_pars.image_size = [m, n];
mask_pars.central_window = [30, 30];
mask_pars.line_num = 70;
mask = load_mask('radial', mask_pars);
% figure;imshow(mask,[]);

% Radial Mask
% mask_pars.image_size = [m, n];
% mask_pars.central_window = [30, 30];
% mask_pars.line_num = 70;
% mask = load_mask('radial', mask_pars);
% figure;imshow(mask,[]);

% Slice Mask
% mask_pars.image_size = [m, n];
% mask_pars.central_window = [30, 30];
% mask_pars.d = 4;
% mask = load_mask('slice', mask_pars);
% figure;imshow(mask,[]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% TVLR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


A = multi_p2DFT(mask, T);
B = A*data;
pars.lambda_1 = lambda_1;
pars.lambda_2 = lambda_2;
pars.max_iter = max_iter;
pars.tol = 1e-10;
pars.verbose = 1;
pars.debug_output = 1;
pars.gt = data;
% for lambda_2 = (2.^[-10:-5])
% pars.lambda_2 = lambda_2;
t0 = tic; [X, psnr_vals, energy_vals, time_vals, rerr, err] = TVLR(A, B, pars); t = toc(t0);
fprintf('RMSE: %.5f, time=%.2fs\n', RMSE(X, data), t);
% end
figure;imshow(abs(sos(X)), []);
figure;plot([1:max_iter+1], err, 'r-');
figure;plot([1:max_iter+1], energy_vals, 'b-');

% A = {}; b = {};
% for t=1:T
% A{t} = p2DFT(mask, size(data(:,:,1)), 1, 2);
% b{t} = A{t}*data(:,:,t);
% im_dc(:,:,t) = A{t}'*b{t};
% end
% input.f=data; input.n1=m;input.n2=n;
% input.alpha=lambda_1;
% input.L=1;
% input.no=50; %%%%  Change no for different iteration numbers
% 
% fprintf('calling the function JTV.....\n');
% t0 = tic;
% out = FISTA_JTV(b, A, input); t = toc(t0);
% im_jtv = out.y;
% figure;imshow(abs(sos(im_jtv)),[]);
% 
% fprintf('RMSE: %.5f, time=%.2fs\n', RMSE(im_jtv, data), t);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Visualization 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure;
% plot(out.rmse);
% legend('JTVMRI');
% xlabel('Iteration');
% ylabel('RMSE');
% 
% figure; hold on;
% subplot(1,5,1); imshow(abs(sos(data)),[]); title('Ground Truth');
% subplot(1,5,2); imshow(abs(sos(im_grappa)),[]); title('GRAPPA');
% subplot(1,5,3); imshow(abs(sos(im_cgspirit)), []); title('CGSPIRiT');
% subplot(1,5,4); imshow(abs(sos(im_cssense)),[]); title('CS\_SENSE');
% subplot(1,5,5); imshow(abs(sos(im_jtv)),[]); title('JTVMRI');


