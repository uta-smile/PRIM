function y=grad2(x)
    [m,n] = size(x);
    y = [zeros(m,1) diff(x,1,2)];
end