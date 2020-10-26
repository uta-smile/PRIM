function out=sparseMRICg(input)

%%% by Nov , 2009
data=input.data; mask=input.mask; 
TVWeight =input. alpha; 	% Weight for TV penalty
xfmWeight = input.beta;	% Weight for Transform L1 penalty
Itnlim = input.num;		% Number of iterations


N = size(data); 	% image Size
DN = size(data); 	% data Size

%generate Fourier sampling operator
FT = p2DFT(input.mask, N, 1, 2);
im_ori = input.im_ori;
% XFM = Wavelet('Daubechies',4,4);	% Wavelet
XFM = input.Phi;
% initialize Parameters for reconstruction
param = init;
param.FT = FT;
param.XFM = XFM;
param.TV = TVOP;
param.data = data;
param.TVWeight =TVWeight;     % TV penalty 
param.xfmWeight = xfmWeight;  % L1 wavelet penalty
param.Itnlim = Itnlim;
param.tt0=cputime;
param.xref=XFM*im_ori;
param.fref=im_ori;
out.im_ori=im_ori;

% % figure(100), imshow(abs(im_dc),[]);drawnow;
% res = XFM*input.im_dc;
res = zeros(N);

rate1=norm(param.xref(:)-res(:), 2)/norm(param.xref(:),2);
% rate2=norm(im_ori(:)-im_dc(:))/norm(im_ori(:),2);

% do iterations
timeTrace=[];errTrace=[];snrTrace=[];
fobjTrace=[];

% timeTrace=[0];errTrace=[1];snrTrace=[0];
% fobjTrace=norm(data, 'fro').^2;

for n=1:input.maxIter
	[res, errTraceTmp, timeTraceTmp, snrTraceTmp, fobjTraceTmp] = fnlCgTrace(res,param);
    errTrace=[errTrace, errTraceTmp];
    timeTrace=[timeTrace, timeTraceTmp];
    snrTrace=[snrTrace, snrTraceTmp];
    fobjTrace=[fobjTrace, fobjTraceTmp];
	im_res = XFM'*res;
% 	figure(100), imshow(abs(im_res),[]), drawnow
end

% out.im_dc=im_dc;
out.im_recover=im_res;
out.errTrace=errTrace;
out.timeTrace=timeTrace;
out.snrTrace=snrTrace;
out.fobjTrace=fobjTrace;
return