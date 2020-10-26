function A = multi_p2DFT( mask, T)
%
%   multi_p2DFT
%
[m, n] = size(mask);
A_single = p2DFT(mask, [m n], 1, 2);
A = A_operator(@(X) multi_A(A_single, X, T, [m, n]), @(X) multi_A(A_single', X, T, [m, n]));

end

function R = multi_A(A_single, X, T, X_single_size)
    R = zeros(size(X));
    N = prod(X_single_size);
    for t = 1:T
        idx = (t-1)*N;
        R(idx+1:idx+N) = A_single*reshape(X(idx+1:idx+N), X_single_size);
    end
end

