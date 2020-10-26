clear all;close all;
load('cardiac_perf_R8.mat');
[nx,ny,nt,nc]=size(kdata);
E = Emat_xyt(kdata(:,:,:,1)~=0,b1);
B = E'*kdata;
