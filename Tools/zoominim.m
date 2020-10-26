function [out, zoomed]= zoominim(image,flag,x,y)

out = [];
[m,n,T] = size(image);


for t=1:T
    F = image(:,:,t);
    
    F2 = F(x:x+31,y:y+31);
    
    if t==1
    F(x:x+31,y) = 255;
    F(x:x+31,y+31) = 255;
    F(x,y:y+31) = 255;
    F(x+31,y:y+31) = 255;
    end
    
    if flag ==1
    F(193:256,1:64) = imresize(F2,[64,64]);
    F(193:256,64) = 255;
    F(193,1:64) = 255;
    end
    
%    out = cat(2,out,F);
    out(:,:,t) = F;
    
    zoomed(:,:,t) = imresize(F2,[64,64]);
end



% figure; hold on;
% subplot(2,3,1); imshow(out(:,:,1),[]); title('Original');
% subplot(2,3,2); imshow(out(:,:,2),[]); title('GRAPPA');
% subplot(2,3,3); imshow(out(:,:,3),[]); title('CGSPIRiT');
% subplot(2,3,4); imshow(out(:,:,4),[]); title('CSSENSE');
% subplot(2,3,5); imshow(out(:,:,5),[]); title('CaLMMRI');
% subplot(2,3,6); imshow(out(:,:,6),[]); title('Proposed');