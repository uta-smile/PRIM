function z = shrink_group_MT(G,x,y,nonneg,groups,g,T)
%
%  Group-wise or row-wise shrinkage operator: Shrink(x,y)
%
if nonneg  % nonnegativity
    x = max(x,0);
end

if nargin == 7 % group-wise
    for t=1:T
    Gx(:,t) = G*x(:,t);
    Gx2(:,t) = Gx(:,t).^2;
    temp(:,t) = g*Gx2(:,t);
    end
    
    
    tp = max(0,1-y./sqrt(sum(temp,2)));

    groups=repmat(groups, [1, T]);
    z = tp(groups).*Gx;   
end


end