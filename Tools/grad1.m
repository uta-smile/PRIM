
function y=grad1(x)
    [m,n] = size(x);
    y = [zeros(1,n);diff(x,1,1)];
end
    
