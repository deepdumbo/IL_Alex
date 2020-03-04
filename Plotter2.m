%{
    last update 2020/03/04, Hossein

    Function to draw plots from QSM images and save them at the specified location.
    Plotter2 is used in Remove_MEDI_influence_runme.m

    Inputs: QSM images
    Output: the following subplots saved in one figures:
        Subplot 1: absolute iField
        Subplot 2: FSL BET Mask Before
        Subplot 5: Mask imposed on iField
        Subplot 6: Negative of mask
        Subplot 7: Unwrap, iFreq, slice #
        Subplot 8: QSM, slice #
        Subplot 9: Bkgd Removal, slice #
        Subplot 10: R2s, slice #
        Subplot 11: CSF Mask, slice
        Subplot 12: QSM with Mask, slice #
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[]=Plotter2(images,FigSaveLocation)


% lims=[-2,2];

figure('Position', [0,0,1920,1080])

%Field
subplot (3,4,1)
imagesc(abs(images.iField(:,:,images.slice)))%,lims)
axis image
title(sprintf('absolute iField, slice %d', images.slice))
colorbar

%Mask
subplot (3,4,2)
imagesc(images.Mask_red(:,:,images.slice))%,lims)
axis image
title(sprintf('FSL BET Mask Before, slice %d', images.slice))
colorbar


%Mask of brain
subplot (3,4,5)
imagesc(images.Impose_Mask(:,:,images.slice))%,lims)
axis image
title(sprintf('Mask imposed on iField, slice %d', images.slice))
colorbar

%Negative mask
subplot (3,4,6)
imagesc(images.Neg_Mask(:,:,images.slice))%,lims)
axis image
title(sprintf('Negative of mask, slice %d', images.slice))
colorbar


%Unwrap
% subplot (3,4,7)
% imagesc(images.iFreq(:,:,images.slice))%,[11,16])
% axis image
% title(sprintf('%s Unwrap, iFreq, slice %d', images.UnwrapType, images.slice))
% colorbar


%Unwrap
subplot (3,4,7)
imagesc((images.iFreq(:,:,images.slice).*images.Mask_red(:,:,images.slice)))%,[11,16])
axis image
title(sprintf('%s Unwrap, iFreq, slice %d', images.UnwrapType, images.slice))
colorbar

RDFmean=mean(mean(mean(images.RDF)));

%Background removal
subplot (3,4,9)
imagesc(images.RDF(:,:,images.slice),[RDFmean-0.3,RDFmean+0.3]);
axis image
title(sprintf('%s, Bkgd Removal, slice %d', images.RDFtype, images.slice))
colorbar

subplot (3,4,10)
imagesc(images.R2s(:,:,images.slice));%,[-10,60])
axis image
title(sprintf('R2s, slice %d', images.slice))
colorbar
% 
subplot (3,4,11)
imagesc(images.Mask_CSF(:,:,images.slice))%,lims)
axis image
title(sprintf('CSF Mask, slice %d', images.slice))
colorbar


if isreal(images.QSM)==0
    images.QSM=real(images.QSM);
    images.PrettyQSM=real(images.PrettyQSM);
    images.QSMtype=sprintf('real %s', images.QSMtype);
end


%QSM

subplot (3,4,8)
imagesc(images.QSM(:,:,images.slice), [-0.3,0.3])
axis image
title(sprintf('%s QSM, slice %d', images.QSMtype, images.slice))
colorbar


subplot (3,4,12)
imagesc(images.PrettyQSM(:,:,images.slice), [-0.3,0.3])
axis image
title(sprintf('%s QSM with Mask, slice %d', images.QSMtype, images.slice))
colorbar

% figure('Position', [0,0,1920,1080])
% imagesc(RDF_QSM_difference(:,:,images.slice))%, [-0.3,0.3])
% axis image
% title(sprintf('RDF(PDF) to QSM (TKD) difference, slice %d', images.slice))
% colorbar

impixelinfo
colormap(bone)

if isfield(images, 'MaskErode') == 1
    saveas(gcf, [FigSaveLocation, images.dataset '_' num2str(images.TU) '_' images.UnwrapType '_' images.RDFtype '_' images.QSMtype])
    print(gcf, [FigSaveLocation, images.dataset '_' num2str(images.TU) '_' images.UnwrapType '_' images.RDFtype '_' images.QSMtype],'-dpng','-r600')
elseif isfield(images, 'MaskDilate') == 1
    saveas(gcf, [FigSaveLocation, images.dataset '_Dilate_' num2str(images.MaskDilate) '_cube_' images.UnwrapType '_' images.RDFtype '_' images.QSMtype])
    print(gcf, [FigSaveLocation, images.dataset '_Dilate_' num2str(images.MaskDilate) '_cube_' images.UnwrapType '_' images.RDFtype '_' images.QSMtype],'-dpng','-r600')
elseif isfield(images, 'FullMask') == 1
    saveas(gcf, [FigSaveLocation, images.dataset '_' 'NoMask_'  images.UnwrapType '_' images.RDFtype '_' images.QSMtype ])
    print(gcf, [FigSaveLocation, images.dataset '_' NoMask_' images.UnwrapType '_' images.RDFtype '_' images.QSMtype],'-dpng','-r600')
else
    sprintf('Could not save plots. Check how the mask is saved?')
end

figure
imagesc(images.PrettyQSM(:,:,images.slice),[-0.3,0.3]);
axis image
title(sprintf('%s QSM', images.QSMtype));%, PDF_TKD.RDFtype, images.slice))
impixelinfo
colorbar
colormap(bone)

if isfield(images, 'MaskErode') == 1
    saveas(gcf, [FigSaveLocation, images.dataset '_' num2str(images.TU) '_' images.UnwrapType '_' images.RDFtype '_' images.QSMtype '_QSM_only'])
    print(gcf, [FigSaveLocation, images.dataset '_' num2str(images.TU) '_' images.UnwrapType '_' images.RDFtype '_' images.QSMtype '_QSM_only'],'-dpng','-r600')
elseif isfield(images, 'MaskDilate') == 1
    saveas(gcf, [FigSaveLocation, images.dataset '_Dilate_' num2str(images.MaskDilate) '_cube_' images.UnwrapType '_' images.RDFtype '_' images.QSMtype '_QSM_only'])
    print(gcf, [FigSaveLocation, images.dataset '_Dilate_' num2str(images.MaskDilate) '_cube_' images.UnwrapType '_' images.RDFtype '_' images.QSMtype '_QSM_only'],'-dpng','-r600')
elseif isfield(images, 'FullMask') == 1
    saveas(gcf, [FigSaveLocation, images.dataset '_' 'NoMask_'  images.UnwrapType '_' images.RDFtype '_' images.QSMtype '_QSM_only'])
    print(gcf, [FigSaveLocation, images.dataset '_' 'NoMask_' images.UnwrapType '_' images.RDFtype '_' images.QSMtype '_QSM_only'],'-dpng','-r600')
else
    sprintf('Could not save plots. Check how the mask is saved?')
end

%     parula     - Blue-green-orange-yellow color map
%     hsv        - Hue-saturation-value color map.x`
%     hot        - Black-red-yellow-white color map.
%     gray       - Linear gray-scale color map.
%     bone       - Gray-scale with tinge of blue color map.
%     copper     - Linear copper-tone color map.
%     pink       - Pastel shades of pink color map.
%     white      - All white color map.
%     flag       - Alternating red, white, blue, and black color map.
%     lines      - Color map with the line colors.
%     colorcube  - Enhanced color-cube color map.
%     vga        - Windows colormap for 16 colors.
%     jet        - Variant of HSV.
%     prism      - Prism color map.
%     cool       - Shades of cyan and magenta color map.
%     autumn     - Shades of red and yellow color map.
%     spring     - Shades of magenta and yellow color map.
%     winter     - Shades of blue and green color map.
%     summer     - Shades of green and yellow color map.
