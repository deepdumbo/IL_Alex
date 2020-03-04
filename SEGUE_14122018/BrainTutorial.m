
addpath('SEGUE');

%% Load data

% load('BrainExample.mat'); %% Load brain phase image and mask
mask=Mask;
% phase=dataScan.phase;
phase=dataScan.phase(:,:,1:size(Mask,3),:);
% phase=phase.*Mask;

figure;
imshow(rot90(squeeze(phase(:,100,:))),[]);
figure;
imshow(rot90(squeeze(mask(:,100,:))),[]);

%% Use SEGUE for phase unwrapping

Inputs.Mask = mask;
Inputs.Phase = phase;

tic; Unwrapped = SEGUE(Inputs); toc

figure;
imshow(rot90(squeeze(Unwrapped(:,100,:))),[-16 16]);