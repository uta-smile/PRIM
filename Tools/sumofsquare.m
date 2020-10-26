function Y = sumofsquare(X,nomflag)

[m,n] = size(X);

X = X.^2;
Y = sum(X,2);
Y = sqrt(Y);
if nomflag
    normY = max(Y) - min(Y);               
%     normY = repmat(normY, [Y 1]);  
    Y = Y./normY;  
    Y=Y*255;
end