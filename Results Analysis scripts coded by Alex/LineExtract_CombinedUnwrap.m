%{
    last update 2020/03/06, Hossein

    This script draws a line on QSM image and plots the phase of the voxels accross that line. 

    input: QSM images at the CombinedUnwrap stage of the pipeline. images of all 5 echos
    output: plots of phase across one line in the images
%}

%% download the data

% close all
% clear
% [file,path] = uigetfile('*.mat');
% if isequal(file,0)
%    disp('User selected Cancel');
%    return
% else
%    load(fullfile(path,file));
% end


FigSaveLocation='D:\AlexEnsworth\CurrentDataAndFigs\Figures\';

PhaseData=images.PhaseData;
% PhaseData=TempUnwrap;
Mask=images.Mask_red;
mask_4D=repmat(Mask, 1,1,1, 5);
PhaseData=mask_4D.*PhaseData;

%% Set ups
set(gcf, 'Position', get(0,'Screensize'))

k=5;
slice=59;
xcoord = ceil((images.matrix_size(1))*0.54);   

x = [xcoord xcoord];
y = [0 images.matrix_size(1)];
limity = [-7 20];

%% Plot the Unwrapped Phase


for i=1:5
    Img=PhaseData(:,:,slice,i);
    % subplot(1,2,1)
    subplot(3,4,2*i-1)
    imagesc(Img);%, [-5, 28]);
    hold on
    line(x,y);
    hold off
    axis image
    title(sprintf('Phase Unwrap using %s, echo %d', images.UnwrapType, i));
    colorbar
    impixelinfo
    colormap(bone)

    % subplot(1,2,2)
    subplot(3,4,2*i)
    improfile(Img,x,y);
    ylabel('Phase (radians)')
    xlim(y)
    ylim(limity)
    title(sprintf('Line profile of unwrapped phase, echo %d', i));
end


%% Plot the combined phase image


PhsImg=images.iFreq(:,:,slice);


% subplot(1,2,1)
subplot(3,4,11)
imagesc(PhsImg);%, [-5, 28]);
hold on
line(x,y);
hold off
axis image
title(sprintf('Combined phase unwrap using %s', images.UnwrapType));
colorbar
impixelinfo
colormap(bone)
    
% subplot(1,2,2)
subplot(3,4,12)
improfile(PhsImg,x,y);
ylabel('Phase (radians)')
xlim(y)
% ylim([-1 14])
ylim(limity)
title(sprintf('Line profile of combined, unwrapped phase'));


%% Save all

%   Save for unwrapped + RDF comparison
% saveas(gcf, [FigSaveLocation, images.dataset '_Unwrapped_line_plot_with_' images.UnwrapType '_all_echoes_and_combined_' images.TU '_slice_' num2str(slice)])
% print(gcf, [FigSaveLocation, images.dataset '_Unwrapped_line_plot_with_' images.UnwrapType '_all_echoes_and_combined_' images.TU '_slice_' num2str(slice)],'-dpng','-r600')

%% Extras


% subplot(1,2,1)
% imagesc(Img,[-4,23]);
% hold on
% line(x2,y2);
% hold off
% axis image
% title(sprintf('Phase Unwrap using %s, Echo number %d', images.UnwrapType, k));%, PDF_TKD.RDFtype, images.slice))
% colorbar
% impixelinfo
% colormap(bone)
%     
% subplot(1,2,2)
% improfile(Img,x2,y2);
% xlim(y2)
% %axis image
% title(sprintf('Line Profile'));%, PDF_TKD.RDFtype, images.slice))

