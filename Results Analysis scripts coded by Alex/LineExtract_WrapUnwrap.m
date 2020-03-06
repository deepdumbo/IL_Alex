%{
    last update 2020/03/06, Hossein

    This script draws a line on QSM image and plots the phase of the voxels accross that line. 

    input: QSM images at wrapunwrap stage of the pipeline. images of all 5 echos
    output: plots of phase across one line in the images
%}

%% download the data
%load('C:\Users\alexe\OneDrive - McGill University\Alex MSc Thesis Research\QSM\Acquired data Dec 13 2019\Saved_MAT_files\run1_Mask_erosion_4_Quality-Guided-Vero_PDF_iSWIMvero.mat')
close all
% clear
% [file,path] = uigetfile('*.mat');
% if isequal(file,0)
%    disp('User selected Cancel');
%    return
% else
%    load(fullfile(path,file));
% end

% images.slice=65;
PhaseData=images.PhaseData;
Mask=images.Mask_red;
mask_4D=repmat(Mask, 1,1,1, 5);
PhaseData=mask_4D.*PhaseData;


%% Set up params
set(gcf, 'Position', get(0,'Screensize'))

k=5;            %Which Echo do you want?
xcoord = ceil((images.matrix_size(1))*0.54);   

x = [xcoord xcoord];
y = [0 images.matrix_size(1)];
%% Plot the raw phase data


% RawImg=images.iFreq_raw;      %Use this if you have the wrapped image saved
%load('RawPhaseNov8.mat')        %Use this if you need to load in (Specific to the data set, be careful)

PreImg=images.iFreq_raw(:,:,images.slice,k);
% subplot(1,2,1)
subplot(2,2,1)
imagesc(PreImg);
hold on
line(x,y);
hold off
axis image
title(sprintf('Wrapped phase using Echo number %d', k));
colorbar
impixelinfo
colormap(bone)
    
% subplot(1,2,2)
subplot(2,2,2)
improfile(PreImg,x,y);
ylabel('Phase (radians)')
xlim(y)
% ylim([-10,12])
title(sprintf('Line profile of wrapped phase'));


% saveas(gcf, ['./Plots/', images.dataset '_Wrapped_line_plot_x_coord_' num2str(xcoord1)])
% print(gcf, ['./Plots/', images.dataset '_Wrapped_line_plot_x_coord_' num2str(xcoord1)],'-dpng','-r600')

%% Plot the Unwrapped Phase

Img=PhaseData(:,:,images.slice,k);

% subplot(1,2,1)
subplot(2,2,3)
imagesc(Img);%, [-15,23]);
hold on
line(x,y);
hold off
axis image
title(sprintf('Phase Unwrap using %s, Echo number %d', images.UnwrapType, k));
colorbar
impixelinfo
colormap(bone)
    
% subplot(1,2,2)
subplot(2,2,4)
improfile(Img,x,y);
ylabel('Phase (radians)')
xlim(y)
% ylim([-15,12])
title(sprintf('Line profile of unwrapped phase'));

% saveas(gcf, ['./Plots/', images.dataset '_' images.UnwrapType '_line_plot_x_coord_' num2str(xcoord1)])
% print(gcf, ['./Plots/', images.dataset '_' images.UnwrapType '_line_plot_x_coord_' num2str(xcoord1)],'-dpng','-r600')

%% Save all

saveas(gcf, ['./Plots/', images.dataset '_Wrapped_line_plot_and_' images.UnwrapType '_unwrap_x_coord_' num2str(xcoord)])
print(gcf, ['./Plots/', images.dataset '_Wrapped_line_plot_and_' images.UnwrapType '_unwrap_x_coord_' num2str(xcoord)],'-dpng','-r600')



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

