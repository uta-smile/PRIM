function y = partialFourier(idx,n,x,mode)
    if mode==1
       % y = P(idx) * FFT(x)
       z = fft(x) / sqrt(n);
       y = z(idx);
    else
       z = zeros(n,1);
       z(idx) = x;
       y = ifft(z) * sqrt(n);
    end
 end % function partialFourier
