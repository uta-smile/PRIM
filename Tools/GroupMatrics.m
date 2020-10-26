function [G,GmatZ,GmatX,groups,Gr,w]=GroupMatrics(GroupIndex)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tree structure in wavelet coefficients (2D).
% input: GroupIndex. ith row are the parents of i.
%        gamma
% output: 
%         G      -- with GX=Z, X is the wavelet coefficients and Z is the
%                   coefficients with nonoverlap groups. 
%         GmatZ  -- sqrt(GmatZ*(Z^2)) is the L21 norm of Z, with group size
%                     2
%         GmatX  --  sqrt(GmatX*(X^2)) is the L21 norm of X, with group
%                     size 2. X is the original coefficients vector.
%         groups --  groups index in Z.  Z(i) belongs to the groups(i)th group.
%Author: xxx, xxx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[m n]= size(GroupIndex);
groupsize = n+1; 

groups = [];
for i=1:groupsize
    groups = [groups;1:m];
end
groups = groups(:);   

if nargin == 2  
    parentIndex = (1:m)';
    parentIndex = parentIndex*groupsize;
end

I = 1:m;
Iy = [I' GroupIndex]; %add self
Iy = Iy';
Iy = Iy(:);  % construct Y coordinates of G

zerorows = find(Iy==0);

if ~isempty(zerorows)
Iy(zerorows)=[]; %delete the parents of 0
groups(zerorows)=[]; %delete the parents of 0 
%parentIndex(zerorows/groupsize) = [];
%parentIndex = parentIndex - length(zerorows);
end


%WeightZ = [groups(parentIndex) J(parentIndex)' ones(rowsG,1)/gamma];
%WeightZ = spconvert(WeightZ);

rowsG = size(Iy,1);
Ix = (1:rowsG)'; % X coordinates of G
G=[Ix Iy ones(rowsG,1)];
G = spconvert(G);

%lenzeros = length(zerorows);
J=1:rowsG; 
GmatZ = [groups J'  ones(rowsG,1)];
GmatZ = spconvert(GmatZ);
%GmatZ(:,parentIndex) = GmatZ(:,parentIndex)/gamma;
%GmatZ = GmatZ.*GmatZ;
GmatX = [groups Iy ones(rowsG,1)]; % P is the parent-children index in X
GmatX = spconvert(GmatX);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Gr = Iy';
gs = sum(GmatX,2);gs=sqrt(gs);
w = ones(3,groups(end));
w(3,:)=gs;
%scaling
sizescaling=length(zerorows);
w(1,1:sizescaling) = 1:sizescaling;w(2,1:sizescaling) = 1:sizescaling;
w(1,(sizescaling+1):end)=(sizescaling+1):2:(length(Gr)-1);
w(2,(sizescaling+1):end)=(sizescaling+2):2:length(Gr);
%WeightX = [groups(parentIndex) Iy(parentIndex) ones(length(parentIndex),1)*(gamma-1)/gamma];
%WeightX = spconvert(WeightX);
%GmatX = GmatX-WeightX;
end