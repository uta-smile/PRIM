function [Im]=sos(iX)

% SOS
%
% Sum of squares of data in x-space
%
% INPUTS:
%	ix:     complex x-space data
%
% PARALLEL MRI TOOLBOX
%
% Santiago Aja-Fernandez, LPI
% www.lpi.tel.uva.es/~santi
% Valladolid, 28/05/2012

Im=sqrt(sum((abs(iX).^2),3));
