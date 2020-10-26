function out = normalize_MRI(x)
% x = abs(x);

[m,n,T] = size(x);

for t=1:T
tp = x(:,:,t);
tpr = abs(tp(:));
% u = mean(tp);
% s = std(tp);
[l, li] = min(tpr);
u = max(tpr);

y = (tp-tp(li))/(u-l);
out(:,:,t) =reshape(y,[m,n]);


end