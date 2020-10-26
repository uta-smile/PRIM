%% Mask Demo
%
%   This demo shows the load_mask function usage
%


%% Parameter Setting 
clear all; clc;
image_size = [256, 256];
central_window = [30, 30];

%% Random Mask

mask_pars.image_size = image_size;
mask_pars.central_window = central_window;
mask_pars.d = 4;
mask = load_mask('random', mask_pars);
figure;imshow(mask,[]);

%% Radial Mask
mask_pars.image_size = image_size;
mask_pars.central_window = central_window;
mask_pars.line_num = 70;
mask = load_mask('radial', mask_pars);
figure;imshow(mask,[]);

%% Slice Mask
mask_pars.image_size = image_size;
mask_pars.central_window = central_window;
mask_pars.d = 4;
mask = load_mask('slice', mask_pars);
figure;imshow(mask,[]);
