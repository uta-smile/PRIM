%% Data Demo
%
%   A simple demo showing the usage of load_db function
%   

%% Initialization
clear all; close all;

%% Raw Brain Data
%
%   Downloaded from http://www.acsu.buffalo.edu/~jlv27/index_files/software_files/NLGRAPPA.zip
%   Introduction Page: http://www.acsu.buffalo.edu/~jlv27/
% 
[data, ~] = load_db(1);
figure; imshow(abs(sos(data)), []);

%% Brain 8CH Data 
%
%   Downloaded From http://www.eecs.berkeley.edu/~mlustig/software/SPIRiT_v0.3.tar.gz
%   Introduction Page: http://www.eecs.berkeley.edu/~mlustig/Software.html
%
db_pars.resize = [256, 256];
[data, ~] = load_db(2, db_pars);
figure; imshow(abs(sos(data)), []);
