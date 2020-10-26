function res = invD(y,imsize)


[res, flag, relres, iter] = lsqr(@afun,y,[],150,[],[],[],imsize);

res = reshape(res,imsize(1),imsize(2));

function res = afun(x, imsize, istranspose)

if nargin==2
    x = reshape(x,imsize(1), imsize(2));
    res = D(x);
else
    res = adjD(x,imsize);
    res = res(:);
end

