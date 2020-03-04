
addpath('SEGUE');

%% Load data

load('HeadAndNeckExample.mat'); %% Load head-and-neck phase image, and water and fat masks

figure;
imshow(rot90(squeeze(phase(:,120,:))),[]);
figure;
imshow(rot90(squeeze(fatmask(:,120,:))),[]);
figure;
imshow(rot90(squeeze(watermask(:,120,:))),[]);

%% Use SEGUE separately in water and fat masks

Inputs.Mask(:,:,:,1) = fatmask;
Inputs.Mask(:,:,:,2) = watermask;
Inputs.Phase = phase;
tic; Unwrapped = SEGUE(Inputs); toc

figure;
imshow(rot90(squeeze(Unwrapped(:,120,:))),[-16 16]);
