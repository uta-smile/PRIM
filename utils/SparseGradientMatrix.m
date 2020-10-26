%{
function [Q1, Q2]=Grad_Mx(N)
Q1 = sparse(N,N);
n=sqrt(N);
for i=1:N
    if rem(i-1,n)==0
        Q1(i,i)=1;
        continue;
    end
    Q1(i,i)=1;
    Q1(i,i-1)=-1;
end

Q2 = sparse(N,N);
n=sqrt(N);
for i=1:N
    if i<=n
        Q2(i,i)=1;
        continue;
    end
    Q2(i,i)=1;
    Q2(i,i-n)=-1;
end
end
%}

function [Q1, Q2]=Grad_Mx(N)

n = sqrt(N);
vec1 = rem(1:N,n);
a = find(vec1==0);
a = a(1:n-1);
a1 = 2:N; a2 = a1-1;
Q1 = sparse(N,N);

[l,r] = size(a1);

% The diag 1
D1 = ones(N,1);
a_d = 1:N;
% The -1
D = -ones(r,1);
tp = [a1',a2',D];

for id = 1:size(a,2)
    num = a(id);
    tp(num,3) = 0;
end
 
tp1 = [a_d',a_d',D1];
tp = [tp;tp1];
Q1 = spconvert(tp);

 
vec_x = n+1:N;
vec_y = 1:N-n;
vec = [vec_x',vec_y',-ones(N-n,1)];
tp2 = [vec;tp1];
Q2 = spconvert(tp2);


end

