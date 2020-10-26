function [L,S] = lps_ist(param)
%
% L+S reconstruction of undersampled dynamic MRI data using iterative
% soft-thresholding of singular values of L and entries of TS
%
% [L,S]=lps_ist(param)
%
% Input variables and reconstruction parameters must be stored in the
% struct param
%
% param.d: undersampled k-t data (nx,ny,nt,nc)
% param.E: data acquisition operator
% param.T: sparsifying transform
% param.lambda_L: nuclear-norm weight
% param.lambda_S: l1-norm weight
%
% Ricardo Otazo (2013)

M=param.E'*param.d;
[nx,ny,nt]=size(M);
M=reshape(M,[nx*ny,nt]);
Lpre=M;S=zeros(nx*ny,nt); 
ite=0;

fprintf('\n ********** L+S reconstruction **********\n')
% iterations
while(1),
	ite=ite+1;
	% low-rank update
	M0=M;
	[Ut,St,Vt]=svd(M-S,0);
	St=diag(SoftThresh(diag(St),St(1)*param.lambda_L));
	L=Ut*St*Vt';
	% sparse update
	S=reshape(param.T'*(SoftThresh(param.T*reshape(M-Lpre,[nx,ny,nt]),param.lambda_S)),[nx*ny,nt]);
	% data consistency
	resk=param.E*reshape(L+S,[nx,ny,nt])-param.d;
	M=L+S-reshape(param.E'*resk,[nx*ny,nt]);
	% L_{k-1} for the next iteration
	Lpre=L;
	% print cost function and solution update
	tmp2=param.T*reshape(S,[nx,ny,nt]);
	fprintf(' ite: %d , cost: %f3, update: %f3\n', ite,norm(resk(:),2)^2+param.lambda_L*sum(diag(St))+param.lambda_S*norm(tmp2(:),1),norm(M(:)-M0(:))/norm(M0(:))); 
	% stopping criteria 
	if (ite > param.nite) || (norm(M(:)-M0(:))<param.tol*norm(M0(:))), break;end
end
L=reshape(L,nx,ny,nt);
S=reshape(S,nx,ny,nt);
end

% soft-thresholding function
function y=SoftThresh(x,p)
y=(abs(x)-p).*x./abs(x).*(abs(x)>p);
y(isnan(y))=0;
end    