%{
    Last update on 2020/03/06, Hossein
    
    This script loads different images, resizes them and montages them all in the same figure. 
    the sections in this script are the following

    1. make the desired figures (or load saved figures)
    2. crop the images by either
        - specififing boundary coordniates
        - viewing the image and drawing the boundary on it
    3. put the images together and save them 
%}

%function[]=Montage_Plotter(images)
% images.slice=65;
close all

%% 1. make the desired figures (or load saved figures)
% figure('Position', [0,0,1920,1080], 'Menu','none','ToolBar','none'); 
% imagesc(OriginalMEDIUnwrap_QSM(:,:,images.slice),[-0.3,0.3])
% axis image
% axis off
% colormap(bone)
% print(gcf, './MontagePlots/image1','-dpng','-r600')
% close
% 
% figure('Position', [0,0,1920,1080], 'Menu','none','ToolBar','none'); 
% imagesc(FixedMEDIUnwrap_QSM(:,:,images.slice),[-0.3,0.3])
% axis image
% axis off
% colormap(bone)
% print(gcf, './MontagePlots/image2','-dpng','-r600')
% close


% figure('Position', [0,0,1920,1080], 'Menu','none','ToolBar','none'); 
% imagesc(SWIMQSM(:,:,images.slice),[-0.3,0.3])%,lims)
% axis image
% axis off
% colormap(bone)
% print(gcf, './MontagePlots/image3','-dpng','-r300')
% close
% 
% figure('Position', [0,0,1920,1080], 'Menu','none','ToolBar','none'); 
% imagesc(TSVDQSM(:,:,images.slice),[-0.3,0.3]);
% axis image
% axis off
% colormap(bone)
% print(gcf, './MontagePlots/image4','-dpng','-r300')
% close
% 
% figure('Position', [0,0,1920,1080], 'Menu','none','ToolBar','none'); 
% imagesc(real(TVSBQSM(:,:,images.slice)),[-0.3,0.3]);
% axis image
% axis off
% colormap(bone)
% print(gcf, './MontagePlots/image5','-dpng','-r300')
% close

%% 2. crop the images by either
% This is the crop that worked for my home comp
% Image1=imcrop(imread('./MontagePlots/image1.png'), [1733 254 2744 2744]);
% Image2=imcrop(imread('./MontagePlots/image2.png'), [1733 254 2744 2744]);
% Image3=imcrop(imread('./MontagePlots/image3.png'), [1733 254 2744 2744]);
% Image4=imcrop(imread('./MontagePlots/image4.png'), [1733 254 2744 2744]);
% Image5=imcrop(imread('./MontagePlots/image5.png'), [1733 254 2744 2744]);

% This is the crop that works for the work comp
% Image1=imcrop(imread('./MontagePlots/image1.png'), [1759 251 2696 2696]);
% Image2=imcrop(imread('./MontagePlots/image2.png'), [1759 251 2696 2696]);
% Image3=imcrop(imread('./MontagePlots/image3.png'), [1759 251 2696 2696]);
% Image4=imcrop(imread('./MontagePlots/image4.png'), [1759 251 2696 2696]);

%This is the crop that works for the work comp on 600dpi
% Image1=imcrop(imread('./MontagePlots/image1.png'), [527 44 908 848]);
% Image2=imcrop(imread('./MontagePlots/image2.png'), [527 44 908 848]);
% Image3=imcrop(imread('./MontagePlots/image3.png'), [527 44 908 848]);

% Home computer, old head data crop
% Image1=imcrop(imread('./MontagePlots/image1.png'), [2336 255 1540 2744]);
% Image2=imcrop(imread('./MontagePlots/image2.png'), [2336 255 1540 2748]);
% Image3=imcrop(imread('./MontagePlots/image3.png'), [2336 255 1540 2748]);

%As an aside - to determine the crop, use this code:
figure('Position', [0,0,1920,1080], 'Menu','none','ToolBar','none'); 
imagesc(OriginalMEDIUnwrap_QSM(:,:,images.slice))
axis image
axis off
colormap(bone)
print(gcf, './MontagePlots/image1','-dpng','-r600')
[J,rect2] = imcrop(imread('./MontagePlots/image1.png'));

%% 3. put the images together and save them
FileSave='QSM_data_comparison';

montage({Image1, Image2, Image3}, 'Size', [1 3]);
print(gcf, ['./MontagePlots/' FileSave], '-dpng','-r600');
saveas(gcf, ['./MontagePlots/' FileSave]); 








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
