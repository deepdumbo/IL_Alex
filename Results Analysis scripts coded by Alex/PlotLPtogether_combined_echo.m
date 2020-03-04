%% load the data

% Nov 8
% SEG=load('C:\Users\alexe\OneDrive - McGill University\Alex MSc Thesis Research\QSM and T1\QSM\Saved_MAT_files\2020-02-21 Looking at unwrapping across different datasets\SEGUE\Nov_8_Mask_erosion_2_SEGUE_PDF_iSWIMvero.mat');
% MLap=load('C:\Users\alexe\OneDrive - McGill University\Alex MSc Thesis Research\QSM and T1\QSM\Saved_MAT_files\2020-02-21 Looking at unwrapping across different datasets\MEDILap\Nov_8_Mask_erosion_2_MEDI Laplacian_PDF_iSWIMvero.mat');
% QGRG=load('C:\Users\alexe\OneDrive - McGill University\Alex MSc Thesis Research\QSM and T1\QSM\Saved_MAT_files\2020-02-21 Looking at unwrapping across different datasets\QGRG\Nov_8_Mask_erosion_2_Quality-Guided-Vero_PDF_iSWIMvero.mat');

% Dec 13 run 1
% SEG=load('C:\Users\alexe\OneDrive - McGill University\Alex MSc Thesis Research\QSM and T1\QSM\Saved_MAT_files\2020-02-21 Looking at unwrapping across different datasets\SEGUE\Dec_13_run1_Mask_erosion_2_SEGUE_PDF_iSWIMvero.mat');
% MLap=load('C:\Users\alexe\OneDrive - McGill University\Alex MSc Thesis Research\QSM and T1\QSM\Saved_MAT_files\2020-02-21 Looking at unwrapping across different datasets\MEDILap\Dec_13_run1_Mask_erosion_2_MEDI Laplacian_PDF_iSWIMvero.mat');
% QGRG=load('C:\Users\alexe\OneDrive - McGill University\Alex MSc Thesis Research\QSM and T1\QSM\Saved_MAT_files\2020-02-21 Looking at unwrapping across different datasets\QGRG\Dec_13_run1_Mask_erosion_2_Quality-Guided-Vero_PDF_iSWIMvero.mat');

% Dec 13 run 2
% SEG=load('C:\Users\alexe\OneDrive - McGill University\Alex MSc Thesis Research\QSM and T1\QSM\Saved_MAT_files\2020-02-21 Looking at unwrapping across different datasets\SEGUE\Dec_13_run2_Mask_erosion_2_SEGUE_PDF_iSWIMvero.mat');
% MLap=load('C:\Users\alexe\OneDrive - McGill University\Alex MSc Thesis Research\QSM and T1\QSM\Saved_MAT_files\2020-02-21 Looking at unwrapping across different datasets\MEDILap\Dec_13_run2_Mask_erosion_2_MEDI Laplacian_PDF_iSWIMvero.mat');
% QGRG=load('C:\Users\alexe\OneDrive - McGill University\Alex MSc Thesis Research\QSM and T1\QSM\Saved_MAT_files\2020-02-21 Looking at unwrapping across different datasets\QGRG\Dec_13_run2_Mask_erosion_2_Quality-Guided-Vero_PDF_iSWIMvero.mat');

% Jan 22
% SEG=load('C:\Users\alexe\OneDrive - McGill University\Alex MSc Thesis Research\QSM and T1\QSM\Saved_MAT_files\2020-02-21 Looking at unwrapping across different datasets\SEGUE\Siemens_Jan22_Mask_erosion_2_SEGUE_PDF_iSWIMvero.mat');
% MLap=load('C:\Users\alexe\OneDrive - McGill University\Alex MSc Thesis Research\QSM and T1\QSM\Saved_MAT_files\2020-02-21 Looking at unwrapping across different datasets\MEDILap\Siemens_Jan22_Mask_erosion_2_MEDI Laplacian_PDF_iSWIMvero.mat');
% QGRG=load('C:\Users\alexe\OneDrive - McGill University\Alex MSc Thesis Research\QSM and T1\QSM\Saved_MAT_files\2020-02-21 Looking at unwrapping across different datasets\QGRG\Siemens_Jan22_Mask_erosion_2_Quality-Guided-Vero_PDF_iSWIMvero.mat');

% Name the file something unique
FileSave='Nov8_combined_line_profile_comparison';

%% Processing

Mask=QGRG.images.Mask_red;
mask_4D=repmat(Mask, 1,1,1, 5);
SEG_iFreq=mask_4D.*SEG.images.iFreq;
MLap_iFreq=mask_4D.*MLap.images.iFreq;
QGRG_iFreq=mask_4D.*QGRG.images.iFreq;


%% Set up params


% k=5;            %Which Echo do you want?
xcoord = ceil((QGRG.images.matrix_size(1))*0.54);   
slice=QGRG.images.slice;

x = [xcoord xcoord];
y = [0 QGRG.images.matrix_size(1)];
limity=[-0.5 3];



UWcomSeg=SEG_iFreq(:,:,slice);
UWcomMlap=MLap_iFreq(:,:,slice);
UWcomQgrg=QGRG_iFreq(:,:,slice);



UWsegLP=improfile(UWcomSeg,x,y);
UWmlapLP=improfile(UWcomMlap,x,y);
UWqgrgLP=improfile(UWcomQgrg,x,y);

%% Plot the Line Profile

figure
set(gcf, 'Position', get(0,'Screensize'))
hold on
% plot(WrappedLP, 'Color','#984ea3', 'linewidth', 4.5);%, 'linestyle', '-.') 
plot(UWsegLP, 'Color','#4daf4a', 'linewidth', 3.5);%, 'linestyle', ':') 
plot(UWmlapLP, 'Color','#377eb8', 'linewidth', 2.5);%, 'linestyle', '--') 
plot(UWqgrgLP, 'Color','#e41a1c', 'linewidth', 1.5) 
legend({'SEGUE','MEDI Laplacian','QG region growing'},'Location','northeast')
hold off
ylabel('Phase (radians)')
axis square
xlim(y)
ylim(limity)
set(gca, 'FontSize', 24)
set(gca, 'XTick',[])
title(sprintf('Combined phase: the %s dataset', QGRG.images.dataset));


%% Save
print(gcf, ['./Plots/' FileSave], '-dpng','-r600');
saveas(gcf, ['./Plots/' FileSave]); 