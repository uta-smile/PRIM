clear all; close all;
rand('state',0); rand('state',0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load rawdata_brain;
im = ifft2c(raw_data);



im = abs(im);

for t =1:8
    y1= grad1(im(:,:,t));
    y2= grad2(im(:,:,t));
    y(:,:,t) = sqrt(y1.^2 + y2.^2);
end

figure;imshow(y(:,:,1),[]);
figure;imshow(y(:,:,2),[]);
figure;imshow(y(:,:,8),[]);
