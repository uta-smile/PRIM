function res = Idct(w,l,NNx,NNy);


Nx = ceil(NNx/l)*l;

Ny = ceil(NNy/l)*l;




res = zeros(Nx,Ny);
count = 0;
for n=1:l:Nx-7
    for m=1:l:Ny-7
        tmp = reshape(w(count*64+1:count*64+64),8,8);
       % w = w(65:end);
        res(n:n+7,m:m+7) = res(n:n+7,m:m+7) +  idct2(tmp);
        count = count+1;
    end
end


res = res(1:NNx, 1:NNy);


