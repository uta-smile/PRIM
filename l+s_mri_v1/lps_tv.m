function [L,S] = lps_tv(param)
%
% L+S reconstruction of undersampled dynamic MRI data using iterative
% soft-thresholding of singular values of L and iterative clipping of
% entries of S (temporal total variation on S)
%
% [L,S]=lps_ist(param)
%
% Input variables and reconstruction parameters must be stored in the
% struct param
%
% param.d: undersampled k-t data (nx,ny,nt,nc)
% param.E: data acquisition operator
% param.lambda_L: nuclear-norm weight
% param.lambda_S: temporal TV weight
%
% Ricardo Otazo (2013)

M=param.E'*param.d;
[nx,ny,nt]=size(M);
M=reshape(M,[nx*ny,nt]);
S=zeros(nx*ny,nt);
z=zeros(nx*ny,nt-1);
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
    % sparse update - tv using clipping
    z=z+0.25*diff(M-L,1,2);
    z=z./abs(z).*reshape(max(min(abs(z(:)),param.lambda_S/2),-param.lambda_S/2),[nx*ny nt-1]);
    z(isnan(z))=0;
    adjDz(:,1)=-z(:,1);
    adjDz(:,2:nt-1)=-diff(z,1,2);
    adjDz(:,nt)=z(:,end);
    S=M-L-adjDz;
    % data consistency
    resk=param.E*reshape(L+S,[nx,ny,nt])-param.d;
    M=L+S-reshape(param.E'*resk,[nx*ny,nt]);
    % print cost function and solution update
    fprintf(' ite: %d , cost: %f3, update: %f3\n', ite,norm(resk(:),2)^2+param.lambda_L*sum(diag(St))+param.lambda_S*norm(diff(S,1,2),1),norm(M(:)-M0(:))/norm(M0(:))); 
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