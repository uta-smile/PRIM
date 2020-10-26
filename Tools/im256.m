function [y]=im256(I)
% this function maps the gray values to the interval 0--255

I=double(I);
y=(I-min(min(I))).*255/(max(max(I))-min(min(I)));
y=round(y);
% y=int8(y);