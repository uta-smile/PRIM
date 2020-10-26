function [PairRate,PathRate] = getForestStat(thetaMap,DecLevel,threshold)

if nargin < 2  % group-wise
    threshold = 0.75;
end

[m,n,C] = size(thetaMap);

if m~=n
    error('Image need to be square');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% wavelet tree structure
[IdxParent, IdxChildren, Ms] = WaveTreeStructure2D(m,DecLevel);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nonzeroparidx = find(IdxParent>0);

%% pair
sumtheta = sum(thetaMap(nonzeroparidx),2)+sum(thetaMap(IdxParent(nonzeroparidx)),2);

Threshold = C*2*threshold;

suc = find(sumtheta>=Threshold|sumtheta<=(C*2-Threshold));    
PairRate = length(suc) / length(sumtheta);


%% path
nonleafidx = IdxChildren(:,1);
allnodes = 1:m*n;
leafindex =  setdiff(allnodes,nonleafidx);
sizeleaf = length(leafindex);
count =zeros(sizeleaf,C);
for k=1:DecLevel
    for i=1:C
        f = thetaMap(:,:,i);
        f = f(leafindex);
        count(:,i) =  count(:,i) + f(:);
    end
    leafindex = IdxParent(leafindex);
end


count = max(count,DecLevel-count);

count = sum(count,2);
succ = find(count>=C*DecLevel*threshold);
PathRate = length(succ)/length(count);