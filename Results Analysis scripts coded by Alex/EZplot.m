% clear
% [file,path] = uigetfile('*.mat');
% if isequal(file,0)
%    disp('User selected Cancel');
%    return
% else
%    load(fullfile(path,file));
% end
close all
%set(gcf, 'Position', get(0,'Screensize'))

% VQGRG = load('C:\Users\alexe\OneDrive - McGill University\Alex MSc Thesis Research\QSM\Acquired data Dec 13 2019\Saved_MAT_files\Nov_8_Mask_erosion_4_Quality-Guided-Vero_PDF_MEDI.mat');
% slice=VQGRG.images.slice;
% Mask=VQGRG.images.Mask_red;
% mask_4D=repmat(Mask, 1,1,1, 5);
% 
% VQGRG_CombPhase=Mask.*VQGRG.images.iFreq;
% 
% Wrapped=mask_4D.*VQGRG.images.iFreq_raw;

% limity=[-2.5 2.5];
% i=1;
% j=6;

% image1 = maps{i,2};
% image2 = maps{j,2};
% dif_image = image1-image2;
% 
% images.RDFtype = 'LBV';
%images.slice=65;

%Vero_phase=phaseModelUW1*delta_TE;

% image1 = QSM;
% image2 = RDF;
% dif_image = image1-image2;

%RDF:PDF
%subplot(2,3,1)
% QGRG=VQGRG_CombPhase(:,:,slice);
% imagesc(images.Residuals(:,:,65), [0 0.001]);
% axis image
% title(sprintf('Combined Phase Residuals for %s', images.UnwrapType));%, PDF_TKD.RDFtype, images.slice))
% impixelinfo
% colorbar
% colormap(gray)
% img=images.Mask_red.*images.Residuals;
% imagesc(img(:,:,images.slice), [0 8]);
% axis image
% title(sprintf('Combined phase residuals for %s', images.UnwrapType));%, PDF_TKD.RDFtype, images.slice))
% impixelinfo
% colorbar
% colormap(gray)

imagesc(images.iMag(:,:,images.slice))%), [-0.2 0.2]);
axis image
% title(sprintf('Combined phase residuals for %s', images.UnwrapType));%, PDF_TKD.RDFtype, images.slice))
impixelinfo
colorbar
colormap(gray)


% subplot(1,3,2)
% imagesc(RDF(:,:,65));%, [-0.3,0.3]);
% axis image
% title(sprintf('RDF'));%, PDF_TKD.RDFtype, images.slice))
% colorbar
% 
% 
% subplot(1,3,3)
% imagesc(QSM(:,:,65), [-0.3,0.3]);
% axis image
% title(sprintf('QSM'));%, PDF_TKD.RDFtype, images.slice))
% colorbar

%colormap(bone)

print(gcf, ['D:\AlexEnsworth\CurrentDataAndFigs\Figures\MagnitudeImage_for_CEGEP_presentation'],'-dpng','-r600')

% print(gcf, ['./figures/Plots/Combination_of_Phase_true_Residuals_for_' images.UnwrapType '_using_' images.dataset],'-dpng','-r600')
% saveas(gcf, ['./figures/Plots/Combination_of_Phase_true_Residuals_for_' images.UnwrapType '_using_' images.dataset])
% saveas(gcf, ['./Plots/Combination_of_Phase_Residuals_for_' images.UnwrapType])

% %QSM:TKD
% subplot(2,4,2)
% imagesc(PDF_TKD.QSM(:,:,images.slice),[-0.3,0.3]);
% axis image
% title(sprintf('%s QSM using %s, slice %d', PDF_TKD.QSMtype, PDF_TKD.RDFtype, images.slice))
% colorbar
% 
% %QSM:MEDI
% subplot(2,4,3)
% imagesc(PDF_MEDI.QSM(:,:,images.slice),[-0.3,0.3]);
% axis image
% title(sprintf('%s QSM using %s, slice %d', PDF_MEDI.QSMtype, PDF_MEDI.RDFtype, images.slice))
% colorbar
% 
% %QSM:SWIM
% subplot(2,4,4)
% imagesc(LBV_SWIM.QSM(:,:,images.slice),[-0.3,0.3]);
% axis image
% title(sprintf('%s QSM using %s, slice %d', PDF_SWIM.QSMtype, PDF_SWIM.RDFtype, images.slice))
% colorbar
% 
% %%%%%
% %RDF:LBV
% subplot(2,4,5)
% imagesc(LBV_TKD.RDF(:,:,images.slice),[-0.07,0.18]);
% axis image
% title(sprintf('RDF using %s, slice %d', LBV_TKD.RDFtype, images.slice))
% colorbar
% 
% %QSM:TKD
% subplot(2,4,6)
% imagesc(LBV_TKD.QSM(:,:,images.slice),[-0.3,0.3]);
% axis image
% title(sprintf('%s QSM using %s, slice %d', LBV_TKD.QSMtype, LBV_TKD.RDFtype, images.slice))
% colorbar
% 
% %QSM:MEDI
% subplot(2,4,7)
% imagesc(LBV_MEDI.QSM(:,:,images.slice),[-0.3,0.3]);
% axis image
% title(sprintf('%s QSM using %s, slice %d', LBV_MEDI.QSMtype, LBV_MEDI.RDFtype, images.slice))
% colorbar
% 
% %QSM:SWIM
% subplot(2,4,8)
% imagesc(LBV_SWIM.QSM(:,:,images.slice),[-0.3,0.3]);
% axis image
% title(sprintf('%s QSM using %s, slice %d', LBV_SWIM.QSMtype, LBV_SWIM.RDFtype, images.slice))
% colorbar
% 
% 
% 
% impixelinfo
% colormap(bone)

