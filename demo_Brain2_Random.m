%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% This script is used to obtain the results in the paper
%%%% "Calibrationless Parallel MRI with Joint Total Variation
%%%% Regularization"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;
rand('state',0); rand('state',0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load brain_8ch;
im = ifft2c(DATA);
I = sos(im);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sigma = 0.01;     % data noise
alpha = 4e-2; beta = 4e-1;  % regularization parameters
d = 4;   %parameter to control sampling ratio  d=3.45 for 33% d=4 for 25% d=4.4 for 20%
maxIter=50; %iteration number for CSSENSE, SPGL1 and JTVMRI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Preparing the data and operators
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear OMEGA A RR b bb;
%%%%%  Prepare for reference image
m=256; n=m;
% [m,n] = size(I);
I=double(I);
I=imresize(I,[m,n]);
f=I(:);
coils = 8;     
T=coils;
NormI = im256(I);
for t=1:T,
    temp = im256(imresize(im(:,:,t),[m,n]));
    F(:,:,t) = abs(temp);
end
[m n, T] = size(F); N = m*n; f=F(:); pn=m*n;
%%%%%  Prepare for wavelet
DecLevel = 4;
wav = daubcqf(2);
WT = @(x) reshape(mdwt(reshape(x,m,n),wav,DecLevel),N,1);
W = @(x) reshape(midwt(reshape(x,m,n),wav,DecLevel),N,1);
Phi= A_operator(@(x) WT(x), @(x) W(x));    % notice Phi=WT

W2 = @(x) midwt(x,wav);WT2 = @(x) mdwt(x,wav);
Phi2= A_operator(@(x) WT2(x), @(x) W2(x));    % notice Phi=WT
%%%%%  Prepare for sampling mask
samplingratio = 0;
for t=1:1,
    [OMEGA{t}] = RandMask_rect(double(m/d),double(n/d),m,n);
    [mask{t}] = RandMask_InverseTransfer(OMEGA{t},m,n);
    mask{t}(m/2-14:m/2+15,n/2-14:n/2+15) = ones(30,30);
    OMEGA{t} = find(fftshift(mask{t})==1);
    k(t) = 2*length(OMEGA{t})+1;
    samplingratio = samplingratio + length(OMEGA{t});
    
end

for t=2:T
    [OMEGA{t}] = [OMEGA{1}];
    [mask{t}] = [mask{1}];
    k(t) = k(1);
end

samplingratio = samplingratio/(N);
figure;imshow(mask{1},[]);


for t=1:T,
    start=1+sum(k(1:t-1)); stop=sum(k(1:t));
    sz(t,:)=[start, stop, pn];
end
%%%%%  Prepare for SPGL1
for t=1:T,
    f2=F(:,:, t); x0((t-1)*pn+1:t*pn,1)=WT(f2(:));
    f0((t-1)*pn+1:t*pn,1)=f2(:);
    R{t} = @(x) A_fhp_rect(x, OMEGA{t}, m, n);
    RT{t} = @(x) At_fhp_rect(x, OMEGA{t}, m, n);
    AO{t} = @(x) R{t}(W(x)); AOT{t} = @(x) WT(RT{t}(x));
    A{t} = A_operator(@(x) AO{t}(x), @(x) AOT{t}(x));
    sigma = 0.01; noise{t} = sigma*randn(k(t),1);
    start=sz(t,1); stop=sz(t,2);
    b(start:stop,1) = A{t}* x0((t-1)*pn+1:t*pn,1)+ noise{t};
end

%%%%%  Prepare for JTV and CSSENSE
for t=1:T,
    f2=F(:,:, t); ff0(:,t)=f2(:); xs(:,t)=WT(ff0(:,t));
    RR{t} = A_operator(@(x) R{t}(x), @(x) RT{t}(x));
    bb{t} = RR{t}*ff0(:,t) + noise{t};
    
    data = zeros(m,n);data(1,1) = bb{t}(1); KK = length(bb{t});
    data(OMEGA{t}) = bb{t}(2:(KK+1)/2) + i*bb{t}((KK+3)/2:KK);
    im_poc{t} = real(ifft2(data))*sqrt(m*n);
end

%%%%%  Prepare for GRAPPA and SPIRiT
DATA = fft2c(F);
MASK = mask{1};
[CalibSize, dcomp] = getCalibSize(MASK);  % get size of calibration area from mask
DATA = DATA + sigma*randn(size(DATA)) + 1i*sigma*randn(size(DATA)) ; % multiply with sampling matrix
DATA = DATA.*repmat(MASK,[1,1,T]);
kSize = [5,5];
CalibTyk = 0.01;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% GRAPPA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('performing traditional GRAPPA reconstruction');
kCalib = crop(DATA,[CalibSize,T]);
t0=cputime();
res_grappa = GRAPPA(DATA,kCalib,kSize,CalibTyk);
t1=cputime();
im_grappa = ifft2c(res_grappa);
im_grappa_err = im_grappa - F;
im1 = sos(im_grappa);
err1 = sos(im_grappa_err);
im1=im256(im1);

SNR(1) = snr(abs(im_grappa),F);

Time(1) = t1-t0;

fprintf('Iter_time = %.2fsec, snr = %f\n',Time(1), SNR(1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% SPIRiT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('performing SPIRiT')
kCalib = crop(DATA,[CalibSize,T]);
kernel = zeros([kSize,T,T]);

[AtA,] = corrMatrix(kCalib,kSize);
t0=cputime();
for t=1:T
    kernel(:,:,:,t) = calibrate(AtA,kSize,T,t,CalibTyk);
end
GOP = SPIRiT(kernel, 'fft',[m,n]);
nIterCG = 10;
ReconTyk = 1e-5;
[res_cg, RELRES,RESVEC,info] = cgSPIRiT(DATA,GOP,nIterCG,ReconTyk, DATA,F);
t1=cputime();

im_cgspirit = ifft2c(res_cg);
im2 = sos(im_cgspirit);
err2 = sos(im_cgspirit-F);

im2=im256(im2);
SNR(2) = info.snr(end);

Time(2) = t1-t0;
fprintf('Iter_time = %.2fsec, snr = %f \n',Time(2), SNR(2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% CSSENSE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('calling the function CSSENSE.....\n');

for t=1:T,
    FT = p2DFT(mask{t}, [m n], 1, 2);im_ori=F(:,:,t);
    input1.alpha=alpha*2; 	% Weight for TV penalty
    input1.beta=beta*2;	% Weight for Transform L1 penalty
    input1.maxIter=maxIter;input1.num=1;
    input1.mask=mask{t}; input1.im_ori = im_ori;
    noise2=randn(m, n); noise2(find(mask{t}==0))=0;
    input1.data =  FT*(im_ori)+sigma*noise2;
    input1.im_dc=im_poc;
    input1.Phi=Phi2;
    
    out=sparseMRICg(input1); out1=out; %im1=abs(out.im_recover);
    img(:,t)=out.im_recover(:);
    snr0(:,t)=out.snrTrace; time0(:,t)=out.timeTrace;
end
im3=reshape(img, [m, n, T]);
err3 = sos(im3-F);

im3 = sos(im3);
im3=im256(im3);

mtsnr=MT_SNR(snr0, ff0);

SNR(3) = mtsnr(end);
time=sum(time0, 2);
Time(3) = time(end);

fprintf('Iter_time = %.2fsec, snr = %f, samp.ratio = %f \n',Time(3), SNR(3), samplingratio);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% CaLMMRI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('calling the function CaLMMRI.....\n');
opts = spgSetParms('verbosity',0);
opts.Phi=Phi;opts.f0=f0; opts.sz=sz; opts.iterations=maxIter;
nGroups = pn; temp=[1:pn]';temp=repmat(temp, [1, T]);groups=temp(:);

[x,r,g,info] = spg_mt(A,b,groups,0,opts);

for t=1:T
    index = (1:N )+ (t-1)*N;
    X4(:,t) = Phi'*x(index);
end
im4=reshape(X4, [m, n, T]);
err4 = sos(im4-F);

SNR(4)=info.Trace_SNR(end);

im4 = sumofsquare(X4,1);
im4 = reshape(im4,[m,n]);
Time(4) = info.Trace_Time(end);


fprintf('Iter_time = %.2fsec, snr = %f, samp.ratio = %f \n',Time(4), SNR(4), samplingratio);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% JTV
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

input.f=f; input.xs=xs; input.n1=m;input.n2=n;
input.alpha=alpha;input.beta=beta;
input.Phi=Phi;input.L=1;
input.maxitr=5000;input.maxiitr= 5000;
input.stopCriterion=2;
input.tolA=1e-8; input.no=maxIter; %%%%  Change no for different iteration numbers
input.l=-inf; input.u=inf;
input1=input;
input1.alpha=alpha*sqrt(T);input1.beta=beta*sqrt(T);
input1.funv=0;
fprintf('calling the function JTV.....\n');
input1=input;
input1.alpha=alpha*sqrt(T);input1.beta=beta*sqrt(T);
input1.funv=0;

out = JTVMRI(bb, RR, ff0, input1);
X5 = out.y;

im5=reshape(X5, [m, n, T]);
err5 = sos(im5-F);

SNR(5)=out.snr(end);

im5 = sumofsquare(X5,1);

Time(5) = out.xtime(end);

im5 = reshape(im5,[m,n]);

fprintf('Iter_time = %.2fsec, snr = %f, samp.ratio = %f \n',Time(5), SNR(5), samplingratio);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cx = 180;
cy =110;

NormI = reshape(NormI,[m,n]);

ims(:,:,1) = NormI;
ims(:,:,2) = im1;
ims(:,:,3) = im2;
ims(:,:,4) = im3;
ims(:,:,5) = im4;
ims(:,:,6) = im5;
zoomedim = zoominim(ims,1,cx,cy);   % zoom in


figure; hold on;
subplot(2,3,1); imshow(zoomedim(:,:,1),[]); title('Original');
subplot(2,3,2); imshow(zoomedim(:,:,2),[]); title('GRAPPA');
subplot(2,3,3); imshow(zoomedim(:,:,3),[]); title('CGSPIRiT');
subplot(2,3,4); imshow(zoomedim(:,:,4),[]); title('CSSENSE');
subplot(2,3,5); imshow(zoomedim(:,:,5),[]); title('CaLMMRI');
subplot(2,3,6); imshow(zoomedim(:,:,6),[]); title('Proposed');


th1 = 0;th2=30;
figure; hold on;
subplot(2,3,2); imshow(err1,[th1,th2]); title('GRAPPA');
subplot(2,3,3); imshow(err2,[th1,th2]); title('CGSPIRiT');
subplot(2,3,4); imshow(err3,[th1,th2]); title('CSSENSE');
subplot(2,3,5); imshow(err4,[th1,th2]); title('CaLMMRI');
subplot(2,3,6); imshow(err5,[th1,th2]); title('Proposed');