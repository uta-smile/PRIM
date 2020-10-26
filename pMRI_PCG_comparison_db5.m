%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is an example for parrallel MRI 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Related papers:
% 
% Chen Chen, Yeqing Li, and Junzhou Huang, "Calibrationless Parallel MRI with Joint Total Variation Regularization", the 
% Annual International Conference on Medical Image Computing and Computer Assisted Intervention (MICCAI), 2013.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Initialization and Parameter Settings
clear all; close all;
db_pars = struct();
db_pars.resize = [256, 256];
db_pars.normalize = 1;
% db_pars.noise = 2e-3;
[data, fdata] = load_db(5, db_pars);
% regularization parameters
% alpha = 1e3; 
alpha = 1e-6;
beta = 4e-1;

maxIter=100;
[m, n, T] = size(fdata); N = m*n;

%% Random Mask Generation
mask_pars.image_size = [m, n];
mask_pars.central_window = [30, 30];
mask_pars.d = 4;
mask = load_mask('random', mask_pars);
figure;imshow(mask,[]);

%% GRAPPA
kSize = [5,5]; CalibTyk = 0.01;
[CalibSize, dcomp] = getCalibSize(mask);  % get size of calibration area from mask
b = fdata.*repmat(mask,[1,1,T]); % sample the original data
fprintf('======== GRAPPA ==========\n');
kCalib = crop(b, [CalibSize,T]);
t0 = tic; res_grappa = GRAPPA(b,kCalib,kSize,CalibTyk); t_grappa=toc(t0);
im_grappa = ifft2c(res_grappa);
figure;imshow(abs(sos(im_grappa)),[]);

fprintf('RMSE: %.8f, SNR: %.6f, time=%.2fs\n', RMSE(im_grappa, data), snr(im_grappa(:), data(:)), t_grappa);

%% cgSPIRiT
fprintf('======== cgSPIRiT ==========\n');
pars.kSize = [5,5];
pars.nIterCG = 30;
pars.CalibTyk = 0.01;  % Tykhonov regularization in the calibration
pars.ReconTyk = 1e-5;
t0 = tic; [im_cgSPIRiT, ~] = runCGSPIRiT(fdata, mask, pars); t_cgSPIRiT = toc(t0);
figure;imshow(abs(sos(im_cgSPIRiT)),[]);

fprintf('RMSE: %.8f, SNR: %.6f, time=%.2fs\n', RMSE(im_cgSPIRiT, data), snr(im_cgSPIRiT(:), data(:)), t_cgSPIRiT);

%% PocsSPIRiT
% fprintf('======== PocsSPIRiT ==========\n');
% pars.kSize = [5,5];
% pars.nIter = 20;
% pars.CalibTyk = 0.01;  % Tykhonov regularization in the calibration
% % pars.wavWeight = 0.0015;
% pars.wavWeight = 0;
% t0 = tic; [im_pocsSPIRiT, ~] = runPocsSPIRiT(fdata, mask, pars); t_pocsSPIRiT = toc(t0);
% figure;imshow(abs(sos(im_pocsSPIRiT)),[]);
% 
% fprintf('RMSE: %.8f, time=%.2fs\n', RMSE(im_pocsSPIRiT, data), t_pocsSPIRiT);

%% CSSENSE
fprintf('======== CSSENSE ==========\n');

im_cssense = zeros([m, n, T]);
t0 = tic;
for t=1:T,
    A = p2DFT(mask, [m n], 1, 2); 
    t_im = data(:, :, t);
    input1.alpha=2e-4; 	% Weight for TV penalty
    input1.beta=0;	% Weight for Transform L1 penalty
    input1.maxIter=20;input1.num=1;
    input1.mask=mask; input1.im_ori = t_im;
    input1.data =  A*t_im;
    input1.im_dc=real(ifft2(t_im))*sqrt(m*n);
    input1.Phi=Wavelet('Daubechies',4, 4); % DecLevel
    
    fprintf('Processing coil %d...\n', t);
    out=sparseMRICg(input1); out1=out; %im1=abs(out.im_recover);
    im_cssense(:, :, t)=out.im_recover;
    % snr0(:,t)=out.snrTrace; time0(:,t)=out.timeTrace;
end
t_cssense = toc(t0);
figure;imshow(abs(sos(im_cssense)),[]);

fprintf('RMSE: %.8f, SNR: %.6f, time=%.2fs\n', RMSE(im_cssense, data), snr(im_cssense(:), data(:)), t_cssense);

%% SAKE-ESPiRIT
fprintf('======== SAKE-ESPiRIT ==========\n');

pars = struct();
pars.ncalib = 48; pars.ksize = [6,6]; % ESPIRiT kernel-window-size 
pars.sakeIter = 100; pars.wnthresh = 1.8;
pars.eigThresh_im = 0.9;
b = fdata .* repmat(mask,[1,1,T]);
t0 = tic; 
im_sake = runSAKE_ESPIRiT(b, pars);
t_sake = toc(t0);

figure;imshow(abs(sos(im_sake)),[]);
fprintf('RMSE: %.8f, SNR: %.6f, time=%.2fs\n', RMSE(im_sake, data), snr(im_sake(:), data(:)), t_sake);

%% FISTA_JTV
A = {}; b = {};
for t=1:T
A{t} = p2DFT(mask, size(data(:,:,1)), 1, 2);
b{t} = A{t}*data(:,:,t);
im_dc(:,:,t) = A{t}'*b{t};
end
input.f=data; input.n1=m;input.n2=n;
input.alpha=1e-4;
input.L=1;
input.no=maxIter; %%%%  Change no for different iteration numbers

fprintf('======== FISTA_JTV ==========\n');
t0 = tic;
out = FISTA_JTV(b, A, input); t_fista_jtv = toc(t0);
im_fista_jtv = out.y;

figure;imshow(abs(sos(im_fista_jtv)),[]);
fprintf('RMSE: %.8f, SNR: %.6f, time=%.2fs\n', RMSE(im_fista_jtv, data), snr(im_fista_jtv(:), data(:)), t_fista_jtv);

%% PCG_JTV
FT = p2DFT(mask, [m, n], 1, 2);  % partial FFT for the first frame

% prepare measurements
A = cell([1, T]);
b = cell([1, T]);
ratio = zeros([1, T]);
for i = 1:T
    b{i} = FT*data(:, :, i);
    b{i} = b{i}(:);
    ratio(i) = length(find(b{i}~=0))/(m*n);
    A{i} = A_operator(@(x) FT*x, @(x) FT'*x);
end

[D1, D2] = discretefinitematrices(m, n);
input.D1 = D1;input.D2=D2;
input.lambda = alpha; % regularization parameter
input.n1 = m;input.n2 = n;
input.no = 6; % number of iterations for outer loop, may need to change for different datasets
input.cgiter = 15;  % number of iterations for inner PCG loop, may need to change for different datasets
input.ratio = ratio;
input.f = data(:);
input.A = A;
input.b = b;
input.l = -inf; input.u = inf;
input.tol = 0;

fprintf('======== PCG_JTV ==========\n');

t0 = tic; out = FIRLS_JTV(input); t_pcg_jtv = toc(t0);
im_pcg_jtv = out.y;

figure;imshow(abs(sos(im_pcg_jtv)),[]);
fprintf('RMSE: %.8f, SNR: %.6f, time=%.2fs\n', RMSE(im_pcg_jtv, data), snr(im_pcg_jtv(:), data(:)), t_pcg_jtv);




%% Visualization
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


