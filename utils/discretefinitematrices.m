function [D1, D2] = discretefinitematrices(m,n)
% input: 
% m,n: the size of the image

%output:
% D1, D2: the finite difference matrices

N = m*n;
I1 = [1:N 2:N];
J1 = [1:N 1:N-1];
V1 = [ones(1,N) -1*ones(1,N-1)];

indexNull = 1:m:N;

D1 = sparse(I1,J1,V1,N,N);

D1(indexNull,:) = zeros(length(indexNull),N);


I2 = [m+1:N m+1:N];
J2 = [m+1:N 1:N-m];
V2 = [ones(1,N-m) -1*ones(1,N-m)];

% indexNull = 1:m:N;
D2 = sparse(I2,J2,V2,N,N);


