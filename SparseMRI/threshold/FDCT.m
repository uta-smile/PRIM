function res = Fdct(x,l);

[NNx,NNy] = size(x);
Nx = ceil(NNx/l)*l;
Ny = ceil(NNy/l)*l;

xx = zeros(Nx,Ny);
xx(1:NNx,1:NNy) = x;



res = zeros(length(1:8:Nx-8)*length(1:8:Ny-8)*64,1);
count = 0;
for n=1:l:Nx-7
    for m=1:l:Ny-7
        tmp = dct2(xx(n:n+7,m:m+7));
        res(count*64+1:count*64+64) =  tmp(:);
        count = count+1;

    end

end



