%{
    last update 2020/03/06, Hossein
    
    a script that extract the following fields from QSM images and sets the related parameters.
    Then, make the figures mentioned bellow: 

    Extracted Fields
        Mask, SEG_iFreq, QGRG_iFreq, SEG_RDF, QGRG_RDF, RDFSeg, RDFQgrg, RDFsegLP, RDFqgrgLP

    Parameters Set
        xcoord, y, x, limity

    figures made
    imagesc(RDFSeg, imagesc(RDFQgrg, plot(RDFsegLP, plot(RDFqgrgLP
%}

%% download the data
%load('C:\Users\alexe\OneDrive - McGill University\Alex MSc Thesis Research\QSM\Acquired data Dec 13 2019\Saved_MAT_files\run1_Mask_erosion_4_Quality-Guided-Vero_PDF_iSWIMvero.mat')
% [file,path] = uigetfile('*.mat');
% if isequal(file,0)
%    disp('User selected Cancel');
%    return
% else
%    load(fullfile(path,file));
% end
% close all

% Nov 8
SEG=load('C:\Users\alexe\OneDrive - McGill University\Alex MSc Thesis Research\QSM and T1\QSM\Saved_MAT_files\2020-02-21 Looking at unwrapping across different datasets\SEGUE\Nov_8_Mask_erosion_2_SEGUE_PDF_iSWIMvero.mat');
% MLap=load('C:\Users\alexe\OneDrive - McGill University\Alex MSc Thesis Research\QSM and T1\QSM\Saved_MAT_files\2020-02-21 Looking at unwrapping across different datasets\MEDILap\Nov_8_Mask_erosion_2_MEDI Laplacian_PDF_iSWIMvero.mat');
QGRG=load('C:\Users\alexe\OneDrive - McGill University\Alex MSc Thesis Research\QSM and T1\QSM\Saved_MAT_files\2020-02-21 Looking at unwrapping across different datasets\QGRG\Nov_8_Mask_erosion_2_Quality-Guided-Vero_PDF_iSWIMvero.mat');

% Dec 13 run 1
% SEG=load('C:\Users\alexe\OneDrive - McGill University\Alex MSc Thesis Research\QSM and T1\QSM\Saved_MAT_files\2020-02-21 Looking at unwrapping across different datasets\SEGUE\Dec_13_run1_Mask_erosion_2_SEGUE_PDF_iSWIMvero.mat');
% % MLap=load('C:\Users\alexe\OneDrive - McGill University\Alex MSc Thesis Research\QSM and T1\QSM\Saved_MAT_files\2020-02-21 Looking at unwrapping across different datasets\MEDILap\Dec_13_run1_Mask_erosion_2_MEDI Laplacian_PDF_iSWIMvero.mat');
% QGRG=load('C:\Users\alexe\OneDrive - McGill University\Alex MSc Thesis Research\QSM and T1\QSM\Saved_MAT_files\2020-02-21 Looking at unwrapping across different datasets\QGRG\Dec_13_run1_Mask_erosion_2_Quality-Guided-Vero_PDF_iSWIMvero.mat');

% Dec 13 run 2
% SEG=load('C:\Users\alexe\OneDrive - McGill University\Alex MSc Thesis Research\QSM and T1\QSM\Saved_MAT_files\2020-02-21 Looking at unwrapping across different datasets\SEGUE\Dec_13_run2_Mask_erosion_2_SEGUE_PDF_iSWIMvero.mat');
% % MLap=load('C:\Users\alexe\OneDrive - McGill University\Alex MSc Thesis Research\QSM and T1\QSM\Saved_MAT_files\2020-02-21 Looking at unwrapping across different datasets\MEDILap\Dec_13_run2_Mask_erosion_2_MEDI Laplacian_PDF_iSWIMvero.mat');
% QGRG=load('C:\Users\alexe\OneDrive - McGill University\Alex MSc Thesis Research\QSM and T1\QSM\Saved_MAT_files\2020-02-21 Looking at unwrapping across different datasets\QGRG\Dec_13_run2_Mask_erosion_2_Quality-Guided-Vero_PDF_iSWIMvero.mat');

% Jan 22
% SEG=load('C:\Users\alexe\OneDrive - McGill University\Alex MSc Thesis Research\QSM and T1\QSM\Saved_MAT_files\2020-02-21 Looking at unwrapping across different datasets\SEGUE\Siemens_Jan22_Mask_erosion_2_SEGUE_PDF_iSWIMvero.mat');
% % MLap=load('C:\Users\alexe\OneDrive - McGill University\Alex MSc Thesis Research\QSM and T1\QSM\Saved_MAT_files\2020-02-21 Looking at unwrapping across different datasets\MEDILap\Siemens_Jan22_Mask_erosion_2_MEDI Laplacian_PDF_iSWIMvero.mat');
% QGRG=load('C:\Users\alexe\OneDrive - McGill University\Alex MSc Thesis Research\QSM and T1\QSM\Saved_MAT_files\2020-02-21 Looking at unwrapping across different datasets\QGRG\Siemens_Jan22_Mask_erosion_2_Quality-Guided-Vero_PDF_iSWIMvero.mat');

% Name the file something unique
FileSave='Nov8_RDF_compare';



Mask=QGRG.images.Mask_red;
SEG_iFreq=Mask.*SEG.images.iFreq;
QGRG_iFreq=Mask.*QGRG.images.iFreq;
SEG_RDF=Mask.*SEG.images.RDF;
QGRG_RDF=Mask.*QGRG.images.RDF;


%% Set up params
slice=QGRG.images.slice;
xcoord = ceil((QGRG.images.matrix_size(1))*0.54);
% xcoord = 168;   %midline is 168
y = [0 QGRG.images.matrix_size(1)];

x = [xcoord xcoord];
limity = [-0.3 0.6];
%% Plot the combined phase image

RDFSeg=SEG_RDF(:,:,slice);
RDFQgrg=QGRG_RDF(:,:,slice);

RDFsegLP=improfile(RDFSeg,x,y);
RDFqgrgLP=improfile(RDFQgrg,x,y);
% RDFdif=RDFsegLP-RDFqgrgLP;







% saveas(gcf, ['./Plots/', images.dataset '_' images.UnwrapType '_line_plot_x_coord_' num2str(xcoord1)])
% print(gcf, ['./Plots/', images.dataset '_' images.UnwrapType '_line_plot_x_coord_' num2str(xcoord1)],'-dpng','-r600')

%% Plot the RDF



set(gcf, 'Position', get(0,'Screensize'))

subplot(1,3,1)
imagesc(RDFSeg, [-1 1]);
hold on
line(x,y);
hold off
axis image
title(sprintf('RDF by %s and %s unwrap', SEG.images.RDFtype, SEG.images.UnwrapType));
colorbar
impixelinfo
colormap(gray)

subplot(1,3,2)
imagesc(RDFQgrg, [-1 1]);
hold on
line(x,y);
hold off
axis image
title(sprintf('RDF by %s and %s unwrap', QGRG.images.RDFtype, QGRG.images.UnwrapType));
colorbar
impixelinfo
colormap(gray)

subplot(1,3,3)
hold on
% plot(WrappedLP, 'Color','#984ea3', 'linewidth', 4.5);%, 'linestyle', '-.') 
plot(RDFsegLP, 'Color','#4daf4a', 'linewidth', 1.5);%, 'linestyle', ':') 
plot(RDFqgrgLP, 'Color','#e41a1c', 'linewidth', 1.5) 
% plot(RDFdif, 'Color','#377eb8', 'linewidth', 1);%, 'linestyle', '--') 
% legend({'Wrapped','SEGUE','QG region growing','Difference'},'Location','northeast')
legend({'SEGUE','QG region growing'},'Location','northeast')
hold off
ylabel('Phase (radians)')
axis square
xlim(y)
ylim(limity)
% set(gca, 'FontSize', 24)
set(gca, 'XTick',[])
% title(sprintf('Echo %d comparison from the %s dataset', k, QGRG.images.dataset));
title(sprintf('RDF comparison from the %s dataset', QGRG.images.dataset));


% saveas(gcf, ['./Plots/', images.dataset '_Wrapped_line_plot_x_coord_' num2str(xcoord1)])
% print(gcf, ['./Plots/', images.dataset '_Wrapped_line_plot_x_coord_' num2str(xcoord1)],'-dpng','-r600')

%% Save all

%   Save for unwrapped + RDF comparison
print(gcf, ['./Plots/' FileSave], '-dpng','-r600');
saveas(gcf, ['./Plots/' FileSave]); 
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

