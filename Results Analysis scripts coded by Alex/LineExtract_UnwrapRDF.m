%% download the data
%load('C:\Users\alexe\OneDrive - McGill University\Alex MSc Thesis Research\QSM\Acquired data Dec 13 2019\Saved_MAT_files\run1_Mask_erosion_4_Quality-Guided-Vero_PDF_iSWIMvero.mat')
[file,path] = uigetfile('*.mat');
if isequal(file,0)
   disp('User selected Cancel');
   return
else
   load(fullfile(path,file));
end
close all
PhaseData=images.PhaseData;
Mask=images.Mask_red;
mask_4D=repmat(Mask, 1,1,1, 5);
PhaseData=mask_4D.*PhaseData;
%% Set up params
set(gcf, 'Position', get(0,'Screensize'))

xcoord = (size(PhaseData,2))/2;
% xcoord = 168;   %midline is 168
ycoord = size(PhaseData,2);

x = [xcoord xcoord];
y = [0 ycoord];
limity = [-0.7 0.7];
%% Plot the combined phase image

PhsImg=images.iFreq(:,:,45);


% subplot(1,2,1)
subplot(2,2,1)
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
subplot(2,2,2)
improfile(PhsImg,x,y);
ylabel('Phase (radians)')
xlim(y)
ylim([-5 30])
% ylim(limity)
title(sprintf('Line profile of unwrapped phase'));

% saveas(gcf, ['./Plots/', images.dataset '_' images.UnwrapType '_line_plot_x_coord_' num2str(xcoord1)])
% print(gcf, ['./Plots/', images.dataset '_' images.UnwrapType '_line_plot_x_coord_' num2str(xcoord1)],'-dpng','-r600')

%% Plot the RDF

RDFImg=images.RDF(:,:,45);

% subplot(1,2,1)
subplot(2,2,3)
imagesc(RDFImg, [-1 1]);
hold on
line(x,y);
hold off
axis image
title(sprintf('RDF using %s bkgd removal', images.RDFtype));
colorbar
impixelinfo
colormap(bone)
    
% subplot(1,2,2)
subplot(2,2,4)
improfile(RDFImg,x,y);
% ylabel('Phase (radians)')
xlim(y)
ylim(limity)
title(sprintf('Line profile of RDF image'));


% saveas(gcf, ['./Plots/', images.dataset '_Wrapped_line_plot_x_coord_' num2str(xcoord1)])
% print(gcf, ['./Plots/', images.dataset '_Wrapped_line_plot_x_coord_' num2str(xcoord1)],'-dpng','-r600')

%% Save all

%   Save for unwrapped + RDF comparison
saveas(gcf, ['./Plots/', images.dataset '_Unwrapped_line_plot_with_' images.UnwrapType '_and_bkd_removed_RDF_with_' images.RDFtype])
print(gcf, ['./Plots/', images.dataset '_Unwrapped_line_plot_with_' images.UnwrapType '_and_bkd_removed_RDF_with_' images.RDFtype],'-dpng','-r600')

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

