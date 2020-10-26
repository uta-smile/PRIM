
clear all; close all;
db_pars.resize = [64, 64];
db_pars.normalize = 1;
[data, fdata] = load_db(1, db_pars);
% [data, fdata] = load_db(2);

[m, n, T] = size(data);

mask_pars.image_size = [m, n];
mask_pars.central_window = [30, 30];
% mask_pars.d = 4;
mask_pars.line_num = 70;
mask = load_mask('radial', mask_pars);


% parameter
max_iter_fista = 100;
lambda_pcg = 1e-5;
lambda_fista = 1e-5;


%% FISTA\_JTV

%inner loop iteration = 1

A = {}; b = {};
for t=1:T
A{t} = p2DFT(mask, size(data(:,:,1)), 1, 2);
b{t} = A{t}*data(:,:,t);
im_dc(:,:,t) = A{t}'*b{t};
end
input.f=data; 
input.n1=m;input.n2=n;
input.alpha=lambda_fista;
input.L=1;
input.no=max_iter_fista; %%%%  Change no for different iteration numbers

fprintf('calling the function JTV.....\n');
t0 = tic;
out = FISTA_JTV(b, A, input); t1 = toc(t0);
im_fista_jtv = out.y;
rmse_fista = out.rmse;
funval_fista = out.funval;
snr_fista = out.snr;
xtime_fista = out.xtime;

fprintf('RMSE: %.8f, time=%.2fs\n', RMSE(im_fista_jtv, data), t1);


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


% initialization

[D1, D2] = Grad_Mx_revised(m*n);
input.D1 = D1;input.D2=D2;
input.lambda = lambda_pcg; % regularization parameter
input.n1 = m;input.n2 = n;
input.no = max_iter_fista; % number of iterations for outer loop, may need to change for different datasets
input.cgiter = 15;  % number of iterations for inner PCG loop, may need to change for different datasets
input.ratio = ratio;
input.f = data(:);
input.A = A;
% input.At = A';
input.b = b;
input.l = -inf; input.u = inf;
input.tol = 0;

t0 = tic;
out = FIRLS_JTV(input);
t2 = toc(t0);
im_firls_jtv = reshape(out.y, [m, n, T]);
rmse_firls = out.rmse;
funval_firls = out.funval;
snr_firls = out.snr;
xtime_firls = out.xtime;

fprintf('RMSE: %.8f, time=%.2fs\n', RMSE(im_firls_jtv, data), t2);

figure; imshow(abs(sos(im_fista_jtv)), []);
figure; imshow(abs(sos(im_firls_jtv)), []);


lw = 3; font_size = 22; xl = [0, 100];

figure; hold on; box on; xlim(xl);
plot([0:length(rmse_fista)-1], rmse_fista, 'b-', 'linewidth', lw);
plot([0:length(rmse_firls)-1], rmse_firls, 'r-', 'linewidth', lw);
legend('FISTA\_JTV', 'PRIM');
xlabel('Iteration');
ylabel('RMSE');
set(gca, 'FontSize', font_size-4);
textobj = findobj('type', 'text');
set(textobj, 'fontsize', font_size);
h_xlabel = get(gca,'XLabel');
set(h_xlabel,'FontSize',font_size); 
h_xlabel = get(gca,'YLabel');
set(h_xlabel,'FontSize',font_size); 

figure; hold on; box on; xlim(xl);
plot([0:length(funval_fista)-1], funval_fista, 'b-', 'linewidth', lw);
plot([0:length(funval_firls)-1], funval_firls, 'r-', 'linewidth', lw);
legend('FISTA\_JTV', 'PRIM');
xlabel('Iteration');
ylabel('Function Value');
set(gca, 'FontSize', font_size-4);
textobj = findobj('type', 'text');
set(textobj, 'fontsize', font_size);
h_xlabel = get(gca,'XLabel');
set(h_xlabel,'FontSize',font_size); 
h_xlabel = get(gca,'YLabel');
set(h_xlabel,'FontSize',font_size); 

figure; hold on; box on; xlim(xl);
plot([0:length(snr_fista)-1], snr_fista, 'b-', 'linewidth', lw);
plot([0:length(snr_firls)-1], snr_firls, 'r-', 'linewidth', lw);
legend('FISTA\_JTV', 'PRIM');
xlabel('Iteration');
ylabel('SNR');
set(gca, 'FontSize', font_size-4);
textobj = findobj('type', 'text');
set(textobj, 'fontsize', font_size);
h_xlabel = get(gca,'XLabel');
set(h_xlabel,'FontSize',font_size); 
h_xlabel = get(gca,'YLabel');
set(h_xlabel,'FontSize',font_size); 


% figure;plot(rmse_jtv);



