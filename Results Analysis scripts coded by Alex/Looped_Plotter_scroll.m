
%{
    last updated on 2020/03/06, Hossein

    This script loads a QSM image, then shows the slide in the image in a way you can scroll through them. 
    It can be used to make gits or multiple plots
%}

%function[]=Looped_Plotter(images)

% Load your maps cell array that you want to plot

% You can uncomment this next bit of code and it will allow you to go
% select the file you wish to plot

%You only need to do this once - once the cell array is in the workspace,
%it will keep using that array



% [file,path] = uigetfile('*.mat');
% if isequal(file,0)
%    disp('User selected Cancel');
%    return
% else
%    load(fullfile(path,file));
% end



% This code can create gifs, or just multiple plots. 
% To make gifs, comment all lines with a 1, uncomment all lines with a 2
% To make multiple plots, comment all lines with a 2, uncomment all lines with a 1

close all

num_loops=size(maps,1);

filename = 'QSMMEDI250to2500.gif';


fig=figure;                                                                                               %2
for i = 1:num_loops
%     figure                                                                                                  %1
    imagesc(maps{i,1}(:,:,80),[-0.3,0.3])
    axis image
    title(sprintf('MEDI QSM, threshold %0.3f',maps{i,2}))
    colorbar
    colormap(bone)
%     print(gcf, ['./Scroll_Plots/GraphCuts_PDF_iSWIM_Threshold_' num2str(1000*maps{i,2}) ],'-dpng','-r300')    %1
    
    frame=getframe(fig);                                                                                  %2 
    im=frame2im(frame);                                                                                   %2

    [A,map]=rgb2ind(im,256);                                                                              %2
    if i == 1                                                                                             %2  
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.4);                                    %2                    
    else                                                                                                  %2
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.4);                               %2
    end                                                                                                   %2
end

close all




% saveas(gcf, ['./Plots/',num2str(images.MaskErode) '_cube_' images.UnwrapType '_' images.RDFtype '_' images.QSMtype '_' images.orient])
% print(gcf, ['./Plots/',num2str(images.MaskErode) '_cube_' images.UnwrapType '_' images.RDFtype '_' images.QSMtype '_' images.orient],'-dpng','-r600')

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
